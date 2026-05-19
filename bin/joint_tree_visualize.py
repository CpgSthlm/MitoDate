#!/usr/bin/env python3
"""
Phylogenetic tree plotting utilities.

Usage
-----
    python joint.tree.plot.py -i joint.tree.out -m joint.meta [-o joint.tree]

Inputs
------
-i / --input      Tree file (Nexus or Newick, e.g. joint.tree.out from BEAST).
-m / --metadata   Joint metadata TSV/CSV with columns:
                      Sample_ID, Species, Origin, Group-By, Age, Age_Uncertainty
                  Colors are assigned per unique value in Group-By.

Outputs
-------
-o / --output     Output stem (default: joint.tree). Produces:
                      <stem>.png, <stem>.pdf
                  (tree + legend in upper-left corner of the same figure)
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from Bio import Phylo


# Color-blind-friendly qualitative palette (10 colors).
# Based on seaborn's "colorblind" palette (an extension of Wong / Okabe-Ito),
# vetted for deuteranopia / protanopia and robust to grayscale printing.
CB_PALETTE: Tuple[str, ...] = (
    "#0173B2",  # blue
    "#DE8F05",  # orange
    "#029E73",  # green
    "#D55E00",  # vermillion
    "#CC78BC",  # pink-purple
    "#CA9161",  # tan
    "#FBAFE4",  # light pink
    "#949494",  # neutral grey
    "#ECE133",  # yellow
    "#56B4E9",  # sky blue
)
NA_COLOR = "#BBBBBB"

EMPTY_GROUP_VALUES = {"", "nan", "NA", "None", "na", "none"}


@dataclass(frozen=True)
class TipFields:
    sample_id: str
    species: str
    origin: str
    group: str
    age: str


class PhylogenyTreeVisualizer:
    def __init__(
        self,
        metadata: Dict[str, str],
        figsize: Tuple[float, float] = (80, 8),
    ):
        self.metadata = metadata
        self.figsize = figsize
        self.group_colors: Dict[str, str] = {}

    @staticmethod
    def parse_tip_fields(tip_name: str) -> Optional[TipFields]:
        cleaned = (tip_name or "").strip("'\"").strip()
        parts = cleaned.rsplit("_", 4)
        if len(parts) != 5:
            return None
        sample_id, species, origin, group, age = parts
        return TipFields(sample_id, species, origin, group, age)

    def infer_group(self, tip_name: str) -> str:
        fields = self.parse_tip_fields(tip_name)
        if fields is None:
            return "NA"
        grp = self.metadata.get(fields.sample_id)
        if grp is None:
            grp = fields.group
        if grp is None or str(grp).strip() in EMPTY_GROUP_VALUES:
            return "NA"
        return str(grp).strip()

    def sample_groups(self, tree) -> Dict[str, str]:
        return {t.name: self.infer_group(t.name or "") for t in tree.get_terminals()}

    TIP_LABEL_FIELDS = ("none", "full", "sample_id", "species", "origin", "group", "age")

    @classmethod
    def _make_tip_label_func(cls, field: str):
        """Return a Phylo.draw-compatible ``label_func`` for the requested field."""
        field = (field or "none").lower()
        if field == "none":
            return lambda _: ""
        if field == "full":
            return lambda c: (c.name or "") if c.is_terminal() else ""

        def _fn(c):
            if not c.is_terminal():
                return ""
            parsed = cls.parse_tip_fields(c.name or "")
            if parsed is None:
                return c.name or ""
            return getattr(parsed, field, "") or ""

        return _fn

    def assign_colors(self, groups: Iterable[str]) -> Dict[str, str]:
        seen: List[str] = []
        for g in groups:
            if g == "NA":
                continue
            if g not in seen:
                seen.append(g)
        seen.sort()
        mapping = {g: CB_PALETTE[i % len(CB_PALETTE)] for i, g in enumerate(seen)}
        mapping["NA"] = NA_COLOR
        self.group_colors = mapping
        return mapping

    @staticmethod
    def print_na_samples(groups: Dict[str, str]) -> None:
        na_samples = [name for name, grp in groups.items() if grp == "NA"]
        if not na_samples:
            return
        print(f"\n=== {len(na_samples)} NA Samples ===")
        for name in na_samples:
            print(f"- {name}")

    @staticmethod
    def sort_tree(tree) -> None:
        tree.ladderize(reverse=True)

        def sort_clade(clade) -> None:
            if clade.is_terminal() or len(clade.clades) <= 1:
                return
            for child in clade.clades:
                sort_clade(child)
            clade.clades.sort(
                key=lambda c: max((t.name for t in c.get_terminals() if t.name), default=""),
                reverse=True,
            )

        sort_clade(tree.root)

    @staticmethod
    def node_positions(tree) -> Dict:
        terminals = list(tree.get_terminals())
        y_pos = {t: i + 1 for i, t in enumerate(terminals)}
        positions: Dict = {}

        def visit(clade):
            x = tree.distance(tree.root, clade)
            if clade.is_terminal():
                y = y_pos[clade]
                positions[clade] = (x, y)
                return (x, y)
            children = [visit(c) for c in clade.clades]
            y = sum(p[1] for p in children) / len(children)
            positions[clade] = (x, y)
            return (x, y)

        visit(tree.root)
        return positions

    @staticmethod
    def posterior_value(node) -> Optional[float]:
        if node.comment:
            match = re.search(r"posterior=([\d.]+)", node.comment)
            if match:
                return float(match.group(1))
        if node.confidence is not None:
            return float(node.confidence)
        return None

    # BEAST node annotations look like: [&height=1.23,height_95%_HPD={0.5,2.0},posterior=0.95]
    _HPD_RE = re.compile(
        r"height_95%_?HPD\s*=\s*[\{\[]\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*[\}\]]"
    )
    _HEIGHT_RE = re.compile(r"\bheight\s*=\s*([-\d.eE+]+)")

    @classmethod
    def hpd_interval(cls, node) -> Optional[Tuple[float, float]]:
        if not node.comment:
            return None
        match = cls._HPD_RE.search(node.comment)
        if not match:
            return None
        try:
            a, b = float(match.group(1)), float(match.group(2))
        except ValueError:
            return None
        return (min(a, b), max(a, b))

    @classmethod
    def node_height(cls, node) -> Optional[float]:
        if not node.comment:
            return None
        match = cls._HEIGHT_RE.search(node.comment)
        if not match:
            return None
        try:
            return float(match.group(1))
        except ValueError:
            return None

    @staticmethod
    def set_tree_linewidth(ax, linewidth: float) -> None:
        for line in ax.lines:
            line.set_linewidth(linewidth)
        for collection in ax.collections:
            if hasattr(collection, "set_linewidth"):
                collection.set_linewidth(linewidth)

    @staticmethod
    def _nice_step(span: float, target_n: int = 6) -> float:
        """Pick a 1/2/2.5/5 x 10^k step that yields roughly ``target_n`` ticks."""
        if span <= 0:
            return 1.0
        raw = span / target_n
        exp = math.floor(math.log10(raw))
        base = 10 ** exp
        for c in (1, 2, 2.5, 5, 10):
            step = c * base
            if span / step <= target_n + 2:
                return step
        return 10 * base

    @classmethod
    def _set_ybp_xaxis(
        cls,
        ax,
        max_depth: float,
        x_min: float = 0.0,
        x_max: Optional[float] = None,
    ) -> None:
        """Flip x-axis labels so right = 0 yBP (present), left = max yBP.

        ``x_min`` lets callers extend the visible range to the left of the
        root (e.g. to accommodate root-HPD bars that reach beyond the MCC
        root age); ``x_max`` extends it to the right of the youngest tip
        (e.g. to give tip dots and tip labels some headroom). Ticks and
        ``xlim`` are both adjusted accordingly.
        """
        if max_depth <= 0:
            return
        left = min(x_min, 0.0)
        right = max_depth if x_max is None else max(x_max, max_depth)
        ybp_max = max_depth - left
        step = cls._nice_step(ybp_max)
        ybp_ticks: List[float] = []
        y = 0.0
        while y <= ybp_max + step * 1e-6:
            ybp_ticks.append(y)
            y += step
        x_positions = [max_depth - y for y in ybp_ticks]

        def _fmt(val: float) -> str:
            if val == 0:
                return "0"
            if val >= 1000 and val % 1000 == 0:
                return f"{val / 1000:g} ka"
            return f"{val:,.0f}"

        ax.set_xticks(x_positions)
        ax.set_xticklabels([_fmt(v) for v in ybp_ticks])
        ax.set_xlim(left, right)
        ax.set_xlabel("Time (yBP)")

    @staticmethod
    def _compute_right_pad(fig, ax, max_depth: float, has_tip_labels: bool) -> float:
        """Return extra data-units to add to the right of ``max_depth``.

        Always leaves a hair of room so tip dots aren't clipped at ``x=max_depth``.
        If tip labels are present, measures their actual rendered width and
        returns enough padding to fit them.
        """
        base_pad = max(max_depth * 0.005, 1e-9)
        if not has_tip_labels:
            return base_pad
        try:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
        except Exception:
            return max(base_pad, max_depth * 0.12)
        inv = ax.transData.inverted()
        label_pad = 0.0
        for txt in ax.texts:
            if not (txt.get_text() or "").strip():
                continue
            bbox_disp = txt.get_window_extent(renderer=renderer)
            bbox_data = bbox_disp.transformed(inv)
            overflow = bbox_data.x1 - max_depth
            if overflow > label_pad:
                label_pad = overflow
        return max(base_pad, label_pad + max_depth * 0.005)

    @staticmethod
    def save_figure(fig, stem: Path, dpi: int, formats: Iterable[str]) -> None:
        stem.parent.mkdir(parents=True, exist_ok=True)
        for fmt in formats:
            out = stem.with_suffix(f".{fmt}")
            fig.savefig(out, dpi=dpi, bbox_inches="tight", transparent=True)
            print(f"Saved: {out}")

    def plot_colored_tips_tree(
        self,
        tree,
        stem: Path,
        tree_linewidth: float = 1.0,
        posterior_threshold: float = 0.8,
        formats: Iterable[str] = ("png", "pdf"),
        legend_order: Optional[List[str]] = None,
        show_hpd: bool = True,
        hpd_color: str = "#4477AA",
        hpd_alpha: float = 0.45,
        hpd_linewidth: float = 4.0,
        tip_labels: str = "none",
        tip_label_fontsize: float = 3.0,
    ) -> None:
        fig, ax = plt.subplots(1, 1, figsize=self.figsize)
        groups = self.sample_groups(tree)
        self.assign_colors(groups.values())
        self.print_na_samples(groups)
        self.sort_tree(tree)
        positions = self.node_positions(tree)

        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            show_confidence=False,
            branch_labels=None,
            label_func=self._make_tip_label_func(tip_labels),
        )
        self.set_tree_linewidth(ax, tree_linewidth)
        if tip_labels != "none":
            # Shift tip labels ~4 pt to the right so the first character doesn't
            # overlap the tip dot (resolution-independent via dpi_scale_trans).
            offset = ScaledTranslation(2 / 72.0, 0.0, fig.dpi_scale_trans)
            for txt in ax.texts:
                txt.set_transform(txt.get_transform() + offset)
                txt.set_fontsize(tip_label_fontsize)

        terminals = list(tree.get_terminals())
        max_depth = max(
            (tree.distance(tree.root, t) for t in terminals), default=0.0
        )

        hpd_drawn = 0
        hpd_x_min = 0.0
        if show_hpd:
            hpd_drawn, hpd_x_min = self._draw_hpd_bars(
                ax, tree, positions,
                max_depth=max_depth,
                color=hpd_color,
                alpha=hpd_alpha,
                linewidth=hpd_linewidth,
            )

        high_nodes = 0
        for node in tree.get_nonterminals():
            post = self.posterior_value(node)
            if post is None or post <= posterior_threshold or node not in positions:
                continue
            x, y = positions[node]
            ax.scatter(x, y, c="black", s=3, zorder=10, marker="o")
            high_nodes += 1

        for i, terminal in enumerate(terminals):
            x = tree.distance(tree.root, terminal)
            grp = groups.get(terminal.name, "NA")
            color = self.group_colors.get(grp, NA_COLOR)
            ax.scatter(x, i + 1, c=color, s=9, zorder=11, edgecolors="black", linewidth=0.5)

        for side in ("top", "left", "right"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_position(("outward", 8))
        ax.tick_params(left=False, labelleft=False)
        ax.set_ylabel("")
        ax.set_title("")

        x_min_pad = hpd_x_min - max_depth * 0.002 if hpd_x_min < 0 else 0.0
        x_max_pad = self._compute_right_pad(
            fig, ax, max_depth=max_depth, has_tip_labels=tip_labels != "none"
        )
        self._set_ybp_xaxis(
            ax, max_depth, x_min=x_min_pad, x_max=max_depth + x_max_pad
        )

        self._add_legend(ax, legend_order=legend_order)

        print("\n=== Posterior Support ===")
        print(f"Nodes with posterior > {posterior_threshold}: {high_nodes}")
        if show_hpd:
            print(f"Nodes with 95% HPD drawn: {hpd_drawn}")

        print("\n=== Sample Distribution ===")
        for grp, count in sorted(Counter(groups.values()).items()):
            print(f"{grp}: {count}")

        self.save_figure(fig, stem=stem, dpi=300, formats=formats)
        plt.close(fig)

    def _draw_hpd_bars(
        self,
        ax,
        tree,
        positions: Dict,
        max_depth: float,
        color: str = "#4477AA",
        alpha: float = 0.45,
        linewidth: float = 4.0,
        zorder: int = 8,
    ) -> Tuple[int, float]:
        """Draw horizontal 95% HPD bars on internal nodes.

        Each node's plotted x position (``x_n``) corresponds to the mean
        (or MCC-annotated) height ``h_n``. Heights are ages relative to the
        most recent tip (BEAST convention). Converting a height ``h`` to an
        x coordinate: ``x = x_n + (h_n - h)``. If the node has no explicit
        ``height=`` annotation we fall back to an ultrametric assumption,
        ``h_n = max_depth - x_n``.
        """
        drawn = 0
        min_x_low = 0.0
        for node in tree.get_nonterminals():
            hpd = self.hpd_interval(node)
            if hpd is None or node not in positions:
                continue
            x_n, y_n = positions[node]
            h_n = self.node_height(node)
            if h_n is None:
                h_n = max_depth - x_n
            h_low, h_high = hpd
            x_low = x_n + (h_n - h_high)
            x_high = x_n + (h_n - h_low)
            ax.hlines(
                y_n, x_low, x_high,
                colors=color, alpha=alpha,
                linewidth=linewidth, zorder=zorder,
            )
            if x_low < min_x_low:
                min_x_low = x_low
            drawn += 1
        return drawn, min_x_low

    def _add_legend(
        self,
        ax,
        legend_order: Optional[List[str]] = None,
    ) -> None:
        if not self.group_colors:
            return

        if legend_order is None:
            legend_order = sorted(g for g in self.group_colors if g != "NA")
            if "NA" in self.group_colors:
                legend_order.append("NA")

        handles = []
        for label in legend_order:
            color = self.group_colors.get(label)
            if color is None:
                continue
            handles.append(
                plt.Line2D(
                    [0], [0],
                    marker="o", color="w",
                    markerfacecolor=color, markersize=10,
                    label=label,
                    markeredgecolor="black", markeredgewidth=1.0,
                )
            )

        if not handles:
            return

        legend = ax.legend(
            handles=handles,
            loc="upper left",
            frameon=True,
            ncol=1,
            fontsize=11,
            handletextpad=0.8,
            borderpad=0.6,
            labelspacing=0.5,
        )
        legend.set_zorder(20)
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("black")
        legend.get_frame().set_linewidth(0.8)


def analyze_tree_structure(tree) -> None:
    print("=== Tree Analysis ===")
    terminals = list(tree.get_terminals())
    max_depth = max(tree.distance(tree.root, t) for t in terminals)
    total_length = tree.total_branch_length()
    print(f"Max depth: {max_depth:.4f}")
    print(f"Terminal nodes: {len(terminals)}")
    print(f"Internal nodes: {len(list(tree.get_nonterminals()))}")
    print(f"Total length: {total_length:.4f}" if total_length else "Total length: N/A")


def load_tree_from_file(filename: str):
    last_error: Optional[Exception] = None
    for fmt in ("nexus", "newick"):
        try:
            tree = Phylo.read(filename, fmt)
            print(f"Loaded {fmt.upper()} tree from {filename}")
            return tree
        except Exception as exc:
            last_error = exc
    print(f"Error loading tree: {last_error}")
    return None


def _sniff_delimiter(path: Path) -> str:
    return "\t" if path.suffix.lower() in (".tsv", ".tab", ".txt", ".meta") else ","


def load_metadata(path: Path) -> Dict[str, str]:
    """Return a mapping Sample_ID -> Group-By value."""
    required = ("Sample_ID", "Group-By")
    primary = _sniff_delimiter(path)
    fallback = "," if primary == "\t" else "\t"
    for delim in (primary, fallback):
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            fields = reader.fieldnames or []
            if all(col in fields for col in required):
                return {
                    row["Sample_ID"].strip(): (row.get("Group-By") or "").strip()
                    for row in reader if row.get("Sample_ID")
                }
    raise ValueError(
        f"Metadata {path} must contain columns {required}."
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot a joint BEAST tree colored by the Group-By metadata column.",
    )
    p.add_argument("-i", "--input", required=True,
                   help="Tree file (Nexus/Newick), e.g. joint.tree.out.")
    p.add_argument("-m", "--metadata", required=True,
                   help="Joint metadata TSV/CSV with Sample_ID and Group-By columns.")
    p.add_argument("-o", "--output", default="joint.tree",
                   help="Output stem (default: joint.tree). "
                        "Produces <stem>.{png,pdf} and <stem>_legend.{png,pdf}.")
    p.add_argument("--figsize", nargs=2, type=float, default=(12, 12),
                   metavar=("W", "H"), help="Tree figure size (default: 12 12).")
    p.add_argument("--posterior-threshold", type=float, default=0.9,
                   help="Posterior support threshold for black node dots (default: 0.9).")
    p.add_argument("--linewidth", type=float, default=1.0,
                   help="Tree branch line width (default: 1.0).")
    p.add_argument("--tip-labels",
                   choices=list(PhylogenyTreeVisualizer.TIP_LABEL_FIELDS),
                   default="sample_id",
                   help="Show tip labels (default: sample_id). 'none' hides "
                        "them; 'full' uses the whole tip name; otherwise show "
                        "just the named field parsed from "
                        "Sample_ID_Species_Origin_Group_Age.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    print("=== Phylogenetic Tree Visualization ===\n")

    tree_file = Path(args.input)
    meta_file = Path(args.metadata)
    stem = Path(args.output)

    metadata = load_metadata(meta_file)
    print(f"Loaded metadata: {len(metadata)} Sample_ID entries from {meta_file}")

    tree = load_tree_from_file(str(tree_file))
    if tree is None:
        return

    analyze_tree_structure(tree)

    viz = PhylogenyTreeVisualizer(metadata=metadata, figsize=tuple(args.figsize))
    print(f"\nGenerating colored tree (legend in upper-left) -> {stem}.{{png,pdf}}")
    viz.plot_colored_tips_tree(
        tree,
        stem=stem,
        tree_linewidth=args.linewidth,
        posterior_threshold=args.posterior_threshold,
        tip_labels=args.tip_labels,
    )
    print("\nPlots saved successfully")


if __name__ == "__main__":
    main()

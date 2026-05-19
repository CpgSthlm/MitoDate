#!/usr/bin/env python3

"""
Generate the joint-tree metadata file, and optionally rename a FASTA file
using that new metadata.

Inputs
------
-i / --input     Collected results CSV from CollectResults.py, with columns:
                     Sample_ID,
                     joint_ESS, prior_ESS, likelihood_ESS,
                     age(root)_mean, age(root)_ESS,
                     clock.rate_mean, clock.rate_ESS,
                     age(Sample)_mean, age(Sample)_Stdev,
                     age(Sample)_HPD_95_Lower, age(Sample)_HPD_95_Upper,
                     age(Sample)_ESS

-m / --metadata  Original metadata table (TSV or CSV) with columns:
                     Sample_ID, Species, Origin, Group-By,
                     Calibrated_yBP, TipDating

-f / --fasta     Input FASTA file. Sequences are renamed by matching
                 Sample_ID and written to --output-fasta.

Outputs
-------
-o / --output          Joint-tree metadata TSV (default: JointTreeMeta.tsv)
                       with columns: Sample_ID, Species, Origin, Group-By,
                                     Age, Age_Uncertainty

-F / --output-fasta    Renamed FASTA (default: joint.renamed.fasta),
                       sequence ids formatted as
                       Sample_ID_Species_Origin_Group-By_Age.

Mapping
-------
    Age             = age(Sample)_mean
    Age_Uncertainty = age(Sample)_Stdev
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


INPUT_MEAN_COL = "age(Sample)_mean"
INPUT_STDEV_COL = "age(Sample)_Stdev"

META_KEEP_COLS = ["Sample_ID", "Species", "Origin", "Group-By"]
OUTPUT_FIELDS = META_KEEP_COLS + ["Age", "Age_Uncertainty"]

EMPTY_GROUP_VALUES = {"", "nan", "NA", "None"}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _to_float(value, context: str = "") -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Cannot convert {value!r} to float ({context}).")


def _sniff_delimiter(path: Path) -> str:
    """Tab for .tsv/.tab/.txt, comma otherwise."""
    return "\t" if path.suffix.lower() in (".tsv", ".tab", ".txt") else ","


def _require_columns(reader: csv.DictReader, required: Iterable[str], path: Path) -> None:
    fieldnames = reader.fieldnames or []
    missing = [c for c in required if c not in fieldnames]
    if missing:
        raise ValueError(
            f"File {path} is missing required columns: {missing}. "
            f"Found: {fieldnames}"
        )


def _format_age_for_name(age) -> str:
    """Format Age for use in a sequence name (no decimals, no spaces)."""
    if age in ("", None):
        return "ND"
    try:
        return str(int(round(float(age))))
    except (TypeError, ValueError):
        return str(age).replace(" ", "_")


# --------------------------------------------------------------------------- #
# Joint-tree metadata generation
# --------------------------------------------------------------------------- #

def load_results(results_path: Path) -> Dict[str, Tuple[float, float]]:
    """Return {sample_id: (age_mean, age_stdev)} from CollectResults output."""
    ages: Dict[str, Tuple[float, float]] = {}
    with results_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter=_sniff_delimiter(results_path))
        _require_columns(reader, (INPUT_MEAN_COL, INPUT_STDEV_COL), results_path)

        for row in reader:
            sample_id = (row.get("Sample_ID") or "").strip()
            if not sample_id:
                continue
            mean = _to_float(row[INPUT_MEAN_COL], context=f"{sample_id} mean")
            stdev = _to_float(row[INPUT_STDEV_COL], context=f"{sample_id} stdev")
            ages[sample_id] = (mean, stdev)
    return ages


def build_rows(
    metadata_path: Path,
    ages: Dict[str, Tuple[float, float]],
) -> Tuple[List[dict], List[str]]:
    """Merge metadata with ages. Return (rows, missing_sample_ids).

    Age resolution per sample:
      1. If the sample appears in `ages` (tip-dated): use age(Sample)_mean / _Stdev.
      2. Else if Calibrated_yBP is a valid number in the metadata: use it with
         Age_Uncertainty = 0.0 (applies to directly-dated and modern samples).
      3. Otherwise: added to `missing_sample_ids` (TipDating=Y sample absent from CSV).
    """
    rows: List[dict] = []
    missing: List[str] = []
    with metadata_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter=_sniff_delimiter(metadata_path))
        _require_columns(reader, META_KEEP_COLS + ["Calibrated_yBP"], metadata_path)

        for row in reader:
            sample_id = (row.get("Sample_ID") or "").strip()
            out = {col: row.get(col, "") for col in META_KEEP_COLS}

            if sample_id and sample_id in ages:
                out["Age"], out["Age_Uncertainty"] = ages[sample_id]
            else:
                calibrated = (row.get("Calibrated_yBP") or "").strip()
                if calibrated and calibrated.upper() != "ND":
                    try:
                        out["Age"] = float(calibrated)
                        out["Age_Uncertainty"] = 0.0
                    except ValueError:
                        out["Age"] = ""
                        out["Age_Uncertainty"] = ""
                        missing.append(sample_id or "<blank>")
                else:
                    out["Age"] = ""
                    out["Age_Uncertainty"] = ""
                    missing.append(sample_id or "<blank>")

            rows.append(out)

    return rows, missing


def write_tsv(output_path: Path, rows: List[dict]) -> None:
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------- #
# FASTA renaming
# --------------------------------------------------------------------------- #

def _is_token_match(haystack: str, needle: str) -> bool:
    """True if `needle` appears in `haystack` bounded by non-alphanumerics."""
    if needle not in haystack:
        return False
    start = haystack.find(needle)
    end = start + len(needle)
    left_ok = (start == 0) or not haystack[start - 1].isalnum()
    right_ok = (end == len(haystack)) or not haystack[end].isalnum()
    return left_ok and right_ok


def find_matching_sample(seq_name: str, sample_ids: Iterable[str]) -> Optional[str]:
    """Find a Sample_ID that matches the sequence name (exact or token match)."""
    sample_ids = list(sample_ids)

    for sid in sample_ids:
        if str(sid) == seq_name:
            return sid

    for sid in sample_ids:
        sid_str = str(sid)
        if _is_token_match(seq_name, sid_str) or _is_token_match(sid_str, seq_name):
            return sid
    return None


def create_new_name(meta_row: dict) -> str:
    """Build new sequence id: Sample_ID_Species_Origin_Group-By_Age."""
    sample_id = str(meta_row.get("Sample_ID", ""))
    species = str(meta_row.get("Species", ""))
    origin = str(meta_row.get("Origin", ""))
    group_by = str(meta_row.get("Group-By", ""))
    if group_by in EMPTY_GROUP_VALUES:
        group_by = "NA"
    age_str = _format_age_for_name(meta_row.get("Age", ""))

    return f"{sample_id}_{species}_{origin}_{group_by}_{age_str}".replace(" ", "_")


def rename_fasta(input_fasta: Path, output_fasta: Path, rows: List[dict]) -> bool:
    """Rename FASTA sequences using the joint-tree metadata rows in memory."""
    try:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
    except ImportError:
        sys.exit(
            "Error: Biopython is required for --fasta. "
            "Install it via `pip install biopython`."
        )

    meta_by_id = {str(r["Sample_ID"]): r for r in rows if r.get("Sample_ID")}
    sample_ids = list(meta_by_id.keys())

    renamed: List = []
    unmatched: List[str] = []
    total = 0

    print(f"Processing FASTA file: {input_fasta}")
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        total += 1
        match = find_matching_sample(record.id, sample_ids)
        if match is None:
            unmatched.append(record.id)
            continue
        new_name = create_new_name(meta_by_id[str(match)])
        renamed.append(SeqRecord(record.seq, id=new_name, description=""))

    if unmatched:
        print(f"\nWarning: {len(unmatched)} sequences not found in metadata")
        for name in unmatched[:5]:
            print(f"  - {name}")
        if len(unmatched) > 5:
            print(f"  ... and {len(unmatched) - 5} more")

    if not renamed:
        print("Error: No sequences were matched and renamed")
        return False

    SeqIO.write(renamed, str(output_fasta), "fasta")
    print(f"\nSuccessfully processed {len(renamed)}/{total} sequences")
    print(f"Renamed FASTA written to: {output_fasta}")

    found = {rec.id.split("_")[0] for rec in renamed}
    missing = set(sample_ids) - found
    if missing:
        print(f"Warning: {len(missing)} metadata samples not found in FASTA")
    else:
        print(f"All {len(sample_ids)} metadata samples were processed")
    return True


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate joint-tree metadata by replacing the last two columns "
            "(Calibrated_yBP, TipDating) of the original metadata with "
            "Age / Age_Uncertainty computed from CollectResults output. "
            "Optionally rename a FASTA file using the new metadata."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Collected results CSV produced by CollectResults.py.",
    )
    parser.add_argument(
        "-m", "--metadata", required=True,
        help="Original metadata file (TSV or CSV) with columns: "
             "Sample_ID, Species, Origin, Group-By, Calibrated_yBP, TipDating.",
    )
    parser.add_argument(
        "-o", "--output", default="JointTreeMeta.tsv",
        help="Output metadata TSV file (default: JointTreeMeta.tsv).",
    )
    parser.add_argument(
        "-f", "--fasta", required=True,
        help="Input FASTA file to be renamed using the new metadata.",
    )
    parser.add_argument(
        "-F", "--output-fasta", default="joint.renamed.fasta",
        help="Output renamed FASTA file (default: joint.renamed.fasta).",
    )
    return parser


def main(argv=None) -> None:
    args = _build_parser().parse_args(argv)

    results_path = Path(args.input)
    metadata_path = Path(args.metadata)
    output_path = Path(args.output)
    fasta_in = Path(args.fasta)
    fasta_out = Path(args.output_fasta)

    for p in (results_path, metadata_path, fasta_in):
        if not p.is_file():
            sys.exit(f"Error: file not found: {p}")

    ages = load_results(results_path)
    rows, missing = build_rows(metadata_path, ages)

    if missing:
        preview = ", ".join(missing[:10])
        more = f" (+{len(missing) - 10} more)" if len(missing) > 10 else ""
        sys.exit(
            f"Error: {len(missing)} sample(s) in {metadata_path.name} have no "
            f"Age in {results_path.name}: {preview}{more}. "
            "Every sample in the joint-tree metadata must have a value."
        )

    write_tsv(output_path, rows)
    print(
        f"Wrote {len(rows)} rows to {output_path} "
        f"(tip-dated ages from {results_path.name}; calibrated ages from {metadata_path.name})."
    )

    out_dir = fasta_out.parent
    if str(out_dir) and not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    if not rename_fasta(fasta_in, fasta_out, rows):
        sys.exit(1)


if __name__ == "__main__":
    main()

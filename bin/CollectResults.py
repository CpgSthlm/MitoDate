#!/usr/bin/env python3

"""
Script to collect results from multiple tip-dating runs and summarize them.

Columns in output CSV:
Sample_ID,
joint_ESS,
prior_ESS,
likelihood_ESS,
age(root)_mean,
age(root)_ESS,
clock.rate_mean,
clock.rate_ESS,
age(Sample_ID)_mean,
age(Sample_ID)_Stdev,
age(Sample_ID)_HPD_95_Lower,
age(Sample_ID)_HPD_95_Upper,
age(Sample_ID)_ESS
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_summary_file(path: Path):
    """
    Parse one summary file and return a dict with all required fields.
    """
    rows = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    # Helper to get a row by exact Parameter name
    def get_row(param_name):
        for r in rows:
            if r.get("Parameter") == param_name:
                return r
        return None

    joint = get_row("joint")
    prior = get_row("prior")
    likelihood = get_row("likelihood")
    age_root = get_row("age(root)")
    clock_rate = get_row("clock.rate")

    # Find the sample age row: Parameter like "age(...)" AND Sample_ID not empty
    sample_age_row = None
    for r in rows:
        param = r.get("Parameter", "")
        sid = r.get("Sample_ID", "")
        if param.startswith("age(") and sid:
            sample_age_row = r
            break

    if sample_age_row is None:
        raise ValueError(
            f"No row with Parameter starting with 'age(' and non-empty Sample_ID found in {path}"
        )

    sample_id = sample_age_row["Sample_ID"]

    def get_float(row, key):
        if row is None:
            return ""
        val = row.get(key, "")
        try:
            return float(val) if val != "" else ""
        except ValueError:
            return val

    result = {
        "Sample_ID": sample_id,
        "joint_ESS": get_float(joint, "ESS"),
        "prior_ESS": get_float(prior, "ESS"),
        "likelihood_ESS": get_float(likelihood, "ESS"),
        "age(root)_mean": get_float(age_root, "Mean"),
        "age(root)_ESS": get_float(age_root, "ESS"),
        "clock.rate_mean": get_float(clock_rate, "Mean"),
        "clock.rate_ESS": get_float(clock_rate, "ESS"),
        "age(Sample)_mean": get_float(sample_age_row, "Mean"),
        "age(Sample)_Stdev": get_float(sample_age_row, "Stdev"),
        "age(Sample)_HPD_95_Lower": get_float(sample_age_row, "HPD_95_Lower"),
        "age(Sample)_HPD_95_Upper": get_float(sample_age_row, "HPD_95_Upper"),
        "age(Sample)_ESS": get_float(sample_age_row, "ESS"),
    }

    return result


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Collect results from multiple tip-dating summary files."
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="Summary CSV files (with header: Parameter,Mean,Stdev,HPD_95_Lower,HPD_95_Upper,ESS,Sample_ID).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="collected_results.csv",
        help="Output CSV file (default: collected_results.csv).",
    )
    args = parser.parse_args(argv)

    input_paths = [Path(p) for p in args.input_files]

    fieldnames = [
        "Sample_ID",
        "joint_ESS",
        "prior_ESS",
        "likelihood_ESS",
        "age(root)_mean",
        "age(root)_ESS",
        "clock.rate_mean",
        "clock.rate_ESS",
        "age(Sample)_mean",
        "age(Sample)_Stdev",
        "age(Sample)_HPD_95_Lower",
        "age(Sample)_HPD_95_Upper",
        "age(Sample)_ESS",
    ]

    results = []
    for p in input_paths:
        try:
            res = parse_summary_file(p)
            results.append(res)
        except Exception as e:
            sys.stderr.write(f"Error parsing {p}: {e}\n")

    with Path(args.output).open("w", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            row["Sample_ID"] = row["Sample_ID"].split("_")[0]
            writer.writerow(row)

    print(f"Wrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()

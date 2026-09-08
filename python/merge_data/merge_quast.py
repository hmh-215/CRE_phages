#!/usr/bin/env python3
"""
merge_quast.py
--------------
Merge QUAST assembly-statistics reports from multiple samples into one
CSV table: one row per sample, one column per QUAST metric (# contigs,
Total length, N50, GC (%), etc). Only the plain-text report.tsv is used
- no PDF/HTML/Icarus visualisation files are touched.

INPUT
-----
A plain text "manifest" file listing one sample per line, in either form:

    /path/to/WS2762512A05.raw.quast
    WS2762512A05<TAB>/path/to/WS2762512A05.raw.quast

The path can point to either:
  - a QUAST output DIRECTORY (the script looks for report.tsv inside it), or
  - a report.tsv FILE directly.

- If only a path is given, the sample name is inferred from the directory
  (or file) name.
- Blank lines and lines starting with '#' are ignored.
- Fields can be separated by a tab or by whitespace.

Each report.tsv is expected in QUAST's standard layout:

    Assembly                     WS2762512A05.contigs.filtered
    # contigs (>= 0 bp)          187
    # contigs (>= 1000 bp)       92
    Total length (>= 0 bp)       5321044
    ...

A report.tsv can also hold MORE THAN ONE assembly column (e.g. if you ran
"quast.py sample1.fasta sample2.fasta ..." together, as in the
"all_samples.raw.quast" combined runs). In that case every assembly column
in the file becomes its own row in the merged output, named
"<manifest_sample>::<assembly_column_name>" so it's traceable back to the
manifest entry it came from.

Metrics differ slightly depending on whether a reference was supplied to
QUAST (e.g. "# misassemblies", "Genome fraction (%)" only appear in
reference-based runs). Samples missing a given metric are simply left
blank in that column - the union of all metrics seen across all inputs is
used to build the header.

USAGE
-----
    python3 merge_quast.py --manifest quast_dirs.txt --output merged_quast.csv

    # keep only a handful of key columns after the merge (comma or space
    # separated, matched by exact QUAST metric name)
    python3 merge_quast.py --manifest quast_dirs.txt --output merged_quast.csv \\
        --keep-metrics "# contigs (>= 0 bp)" "Total length (>= 0 bp)" "N50" "GC (%)"
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def infer_sample_name(path: Path) -> str:
    """Guess a sample name from a directory or report.tsv file path."""
    if path.name.lower() in ("report.tsv", "report.txt"):
        # e.g. .../WS2762512A05.raw.quast/report.tsv -> WS2762512A05.raw.quast
        return path.parent.name
    return path.stem if path.suffix else path.name


def parse_manifest(manifest_path: Path):
    """Return list of (sample_name, path) tuples. path may be a dir or a report.tsv file."""
    entries = []
    with open(manifest_path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) == 1:
                p = Path(parts[0])
                sample = infer_sample_name(p)
            elif len(parts) >= 2:
                sample, p = parts[0], Path(parts[1])
            else:
                print(f"WARNING: manifest line {lineno} could not be parsed, skipping: {line!r}",
                      file=sys.stderr)
                continue
            entries.append((sample, p))
    return entries


def resolve_report_tsv(path: Path) -> Path | None:
    """Given a directory or file path, return the report.tsv path, or None if not found."""
    if path.is_file():
        return path
    if path.is_dir():
        candidate = path / "report.tsv"
        if candidate.exists():
            return candidate
        # fall back to report.txt (older QUAST versions / --output-txt-only)
        candidate_txt = path / "report.txt"
        if candidate_txt.exists():
            return candidate_txt
    return None


def coerce_value(val: str):
    """Try to turn a QUAST report value into a number; keep as string if it isn't one."""
    val = val.strip()
    if val in ("", "-", "NA", "N/A"):
        return None
    try:
        if "." in val:
            return float(val)
        return int(val)
    except ValueError:
        return val


def load_report(sample: str, report_path: Path) -> dict:
    """
    Parse one report.tsv. Returns {row_label: {metric: value}} where
    row_label is the manifest sample name (single-assembly report) or
    "sample::assembly_col" (multi-assembly report).
    """
    with open(report_path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]

    if not lines:
        print(f"WARNING: empty report file for sample '{sample}': {report_path} - skipping",
              file=sys.stderr)
        return {}

    header = lines[0].split("\t")
    if header[0].strip().lower() != "assembly":
        print(f"WARNING: unexpected header in {report_path} (expected first column 'Assembly') "
              f"- skipping sample '{sample}'", file=sys.stderr)
        return {}

    assembly_cols = header[1:]
    if not assembly_cols:
        print(f"WARNING: no assembly columns found in {report_path} - skipping sample '{sample}'",
              file=sys.stderr)
        return {}

    multi = len(assembly_cols) > 1
    row_labels = [f"{sample}::{col}" if multi else sample for col in assembly_cols]

    data = {label: {} for label in row_labels}
    for line in lines[1:]:
        fields = line.split("\t")
        metric = fields[0].strip()
        values = fields[1:]
        # pad in case a row has fewer values than assembly columns
        values += [""] * (len(assembly_cols) - len(values))
        for label, raw_val in zip(row_labels, values):
            data[label][metric] = coerce_value(raw_val)

    return data


def build_merged_table(entries) -> pd.DataFrame:
    all_rows = {}       # row_label -> {metric: value}
    metric_order = []   # preserves first-seen order of metrics
    seen_metrics = set()
    row_order = []

    for sample, path in entries:
        report_path = resolve_report_tsv(path)
        if report_path is None:
            print(f"WARNING: could not find report.tsv/report.txt for sample '{sample}' "
                  f"at {path} - skipping", file=sys.stderr)
            continue

        parsed = load_report(sample, report_path)
        for label, metrics in parsed.items():
            if label in all_rows:
                print(f"WARNING: duplicate row label '{label}' - later entry overwrites earlier one",
                      file=sys.stderr)
            all_rows[label] = metrics
            row_order.append(label)
            for m in metrics:
                if m not in seen_metrics:
                    seen_metrics.add(m)
                    metric_order.append(m)

    # de-duplicate row_order while preserving order
    seen_rows = set()
    ordered_rows = []
    for r in row_order:
        if r not in seen_rows:
            ordered_rows.append(r)
            seen_rows.add(r)

    records = []
    for label in ordered_rows:
        rec = {"Sample": label}
        metrics = all_rows[label]
        for m in metric_order:
            rec[m] = metrics.get(m, None)
        records.append(rec)

    df = pd.DataFrame.from_records(records, columns=["Sample"] + metric_order)
    return df


def main():
    ap = argparse.ArgumentParser(
        description="Merge QUAST report.tsv files from multiple samples into one CSV table."
    )
    ap.add_argument("--manifest", required=True, type=Path,
                     help="Text file listing sample QUAST directories or report.tsv paths "
                          "(one per line, optionally 'sample<TAB>path').")
    ap.add_argument("--output", required=True, type=Path,
                     help="Output CSV path.")
    ap.add_argument("--keep-metrics", nargs="+", default=None,
                     help="Optional: only keep these metric columns (exact QUAST names), "
                          "in the order given. Default: keep every metric encountered.")
    args = ap.parse_args()

    if not args.manifest.exists():
        sys.exit(f"ERROR: manifest file not found: {args.manifest}")

    entries = parse_manifest(args.manifest)
    if not entries:
        sys.exit(f"ERROR: no valid sample entries parsed from {args.manifest}")

    print(f"Parsed {len(entries)} sample entries from manifest.")

    merged = build_merged_table(entries)

    if args.keep_metrics:
        missing = [m for m in args.keep_metrics if m not in merged.columns]
        if missing:
            print(f"WARNING: requested --keep-metrics not found in any report and will be "
                  f"dropped: {missing}", file=sys.stderr)
        keep = ["Sample"] + [m for m in args.keep_metrics if m in merged.columns]
        merged = merged[keep]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)

    n_samples = len(merged)
    n_metrics = merged.shape[1] - 1
    print(f"Merged table: {n_samples} samples x {n_metrics} metrics")
    print(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()

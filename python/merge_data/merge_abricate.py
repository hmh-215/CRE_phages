#!/usr/bin/env python3
"""
merge_abricate.py
------------------
Merge Abricate per-sample TSV outputs into one wide table:

    GENE | DATABASE | ACCESSION | PRODUCT | RESISTANCE | sample1 | sample2 | ...

Each sample column holds the NUMBER OF COPIES of that gene found in that
sample (0 if absent, 1 if found once, 2+ if found more than once - e.g.
duplicated AMR genes, multiple contigs, etc). This mirrors merge_amrfinder.py,
but reads Abricate's column schema instead of AMRFinderPlus's.

INPUT
-----
A plain text "manifest" file listing one sample per line, in either form:

    /path/to/WS2762512A05.abricate.card.tsv
    WS2762512A05<TAB>/path/to/WS2762512A05.abricate.card.tsv

- If only a path is given, the sample name is inferred from the filename
  (common suffixes like ".abricate.card.tsv", "_abricate.tsv", ".tsv" are
  stripped automatically).
- Blank lines and lines starting with '#' are ignored.
- Fields can be separated by a tab or by whitespace.

Each referenced file is expected to be a standard Abricate per-sample output
TSV (the default 15-column format, i.e. the *.abricate.<db>.tsv file, NOT the
`abricate --summary` table), i.e. it must contain at least the columns:
"GENE", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE".

Because PA_comparative-style pipelines run Abricate once per database
(vfdb / card / ncbi / resfinder) and write one TSV per sample per database,
a manifest normally points at a single database's files (e.g.
*.abricate.card.tsv for every sample). --database is provided as a safety
filter in case a manifest accidentally mixes databases.

USAGE
-----
    # Merge CARD results across all samples listed in a manifest
    python3 merge_abricate.py --manifest card_samples.txt --output merged_card.tsv

    # comma-separated output instead of tab-separated
    python3 merge_abricate.py --manifest card_samples.txt --output merged_card.csv --sep ,

    # extra safety net: only keep rows whose DATABASE column says "card"
    # (useful if the manifest was built from a mixed-database glob)
    python3 merge_abricate.py --manifest all_samples.txt --output merged_card.tsv --database card
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Columns Abricate always writes; these are the gene-description columns we
# carry through to the merged table (order preserved). Per-hit-only columns
# (SEQUENCE, START, END, STRAND, COVERAGE, COVERAGE_MAP, GAPS, %COVERAGE,
# %IDENTITY) are intentionally dropped, same as AMRFinderPlus's HMM/identity
# columns are dropped in merge_amrfinder.py - they vary hit-to-hit and don't
# belong in a gene x sample summary.
DESC_COLS = ["GENE", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"]

# Suffixes to strip off a filename when a sample name isn't given explicitly.
NAME_SUFFIXES = [
    ".abricate.vfdb.tsv",
    ".abricate.card.tsv",
    ".abricate.ncbi.tsv",
    ".abricate.resfinder.tsv",
    "_abricate.tsv",
    ".abricate.tsv",
    ".tsv",
    ".txt",
]


def infer_sample_name(path: Path) -> str:
    name = path.name
    for suf in NAME_SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    return path.stem


def parse_manifest(manifest_path: Path):
    """Return list of (sample_name, file_path) tuples."""
    entries = []
    with open(manifest_path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # accept tab or generic whitespace separation
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) == 1:
                fpath = Path(parts[0])
                sample = infer_sample_name(fpath)
            elif len(parts) >= 2:
                sample, fpath = parts[0], Path(parts[1])
            else:
                print(f"WARNING: manifest line {lineno} could not be parsed, skipping: {line!r}",
                      file=sys.stderr)
                continue
            entries.append((sample, fpath))
    return entries


def load_sample_tsv(sample: str, fpath: Path, database_filter: str | None) -> pd.DataFrame:
    if not fpath.exists():
        print(f"WARNING: file not found for sample '{sample}': {fpath} - skipping", file=sys.stderr)
        return pd.DataFrame(columns=DESC_COLS)

    df = pd.read_csv(fpath, sep="\t", dtype=str)
    # Abricate writes its first column as "#FILE" - normalise so callers can
    # rely on the column name without the leading '#'.
    df = df.rename(columns={"#FILE": "FILE"})

    missing = [c for c in DESC_COLS if c not in df.columns]
    if missing:
        print(f"WARNING: sample '{sample}' ({fpath}) is missing expected column(s) "
              f"{missing} - skipping this sample", file=sys.stderr)
        return pd.DataFrame(columns=DESC_COLS)

    if df.empty:
        # Abricate still writes a header-only TSV when no genes are found for
        # a sample; that's valid (0 copies for everything), not a warning.
        return pd.DataFrame(columns=DESC_COLS)

    if database_filter:
        df = df[df["DATABASE"].str.lower() == database_filter.lower()]

    return df[DESC_COLS].copy()


def build_merged_table(entries, database_filter: str | None) -> pd.DataFrame:
    sample_names = []
    per_sample_counts = {}
    desc_lookup = {}  # gene symbol -> description row (first seen wins)
    row_order = []     # preserves first-seen order of genes

    for sample, fpath in entries:
        if sample in per_sample_counts:
            print(f"WARNING: duplicate sample name '{sample}' in manifest - "
                  f"later entry will overwrite the earlier count column", file=sys.stderr)
        sample_names.append(sample)

        df = load_sample_tsv(sample, fpath, database_filter)

        # count copies per gene symbol within this sample
        counts = df.groupby("GENE").size()
        per_sample_counts[sample] = counts

        # record description columns the first time we see each gene symbol
        for _, row in df.drop_duplicates(subset="GENE").iterrows():
            sym = row["GENE"]
            if sym not in desc_lookup:
                desc_lookup[sym] = row
                row_order.append(sym)

    # de-duplicate sample_names while preserving order (last-write-wins handled above)
    seen = set()
    ordered_samples = []
    for s in sample_names:
        if s not in seen:
            ordered_samples.append(s)
            seen.add(s)

    # assemble final table
    records = []
    for sym in row_order:
        rec = {col: desc_lookup[sym][col] for col in DESC_COLS}
        for sample in ordered_samples:
            rec[sample] = int(per_sample_counts[sample].get(sym, 0))
        records.append(rec)

    merged = pd.DataFrame.from_records(records, columns=DESC_COLS + ordered_samples)

    # sort genes alphabetically by symbol for a stable, readable output
    merged = merged.sort_values("GENE", kind="stable").reset_index(drop=True)

    return merged


def main():
    ap = argparse.ArgumentParser(
        description="Merge Abricate per-sample TSVs into one gene x sample copy-count table."
    )
    ap.add_argument("--manifest", required=True, type=Path,
                     help="Text file listing sample paths (one per line, optionally 'sample<TAB>path'). "
                          "Point this at one database's per-sample TSVs at a time "
                          "(e.g. every *.abricate.card.tsv) for a clean merge.")
    ap.add_argument("--output", required=True, type=Path,
                     help="Output table path.")
    ap.add_argument("--sep", default="\t",
                     help="Field separator for the output table (default: tab). Use ',' for CSV.")
    ap.add_argument("--database", default=None, metavar="DB",
                     help="Optional safety filter: only keep rows whose 'DATABASE' column matches "
                          "this value (case-insensitive), e.g. --database card. Default: keep "
                          "whatever database(s) are present in each file as-is. Use this if a "
                          "manifest might mix files from more than one Abricate database "
                          "(vfdb/card/ncbi/resfinder).")
    args = ap.parse_args()

    if not args.manifest.exists():
        sys.exit(f"ERROR: manifest file not found: {args.manifest}")

    entries = parse_manifest(args.manifest)
    if not entries:
        sys.exit(f"ERROR: no valid sample entries parsed from {args.manifest}")

    print(f"Parsed {len(entries)} sample entries from manifest.")
    if args.database:
        print(f"Filtering to DATABASE == '{args.database}'")

    merged = build_merged_table(entries, args.database)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep=args.sep, index=False)

    n_genes = len(merged)
    n_samples = merged.shape[1] - len(DESC_COLS)
    print(f"Merged table: {n_genes} genes x {n_samples} samples")
    print(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()

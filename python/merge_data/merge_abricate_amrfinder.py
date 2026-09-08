#!/usr/bin/env python3
"""
merge_abricate_amrfinder.py
----------------------------
Reconcile a merged Abricate table (from merge_abricate.py) and a merged
AMRFinderPlus table (from merge_amrfinder.py) into one combined gene x sample
table, tagged with a "source" column showing whether each gene was called by
ABRICATE_CARD, AMRFINDER, or both.

BACKGROUND
----------
Abricate and AMRFinderPlus use different reference databases and different
gene-naming conventions for the same underlying resistance gene (e.g. Abricate/
CARD often reports "blaNDM-1" while AMRFinderPlus reports "NDM-1"). A naive
merge on the raw gene-symbol string therefore treats the same gene as two
different genes. This script normalises both sources' gene symbols before
matching, so a supervisor can see at a glance which annotation tool(s)
support a given resistance call - the same true-positive-style reconciliation
idea as prokka_rgi_reconcile.py, but applied across two independent AMR gene
callers instead of Prokka vs RGI.

GENE-NAME MATCHING STRATEGY
----------------------------
1. normalise(): lowercase, strip hyphens/underscores/spaces (same helper as
   prokka_rgi_reconcile.py). "blaNDM-1" -> "blandm1".
2. If --keep-bla-prefix is NOT set (the default), a leading "bla" is then
   stripped from the normalised string before comparison, so "blandm1"
   (Abricate) and "ndm1" (AMRFinder) both resolve to the matching key "ndm1".
   This is a deliberately narrow, single, well-understood rule (the "bla"
   beta-lactamase prefix convention) rather than open-ended fuzzy/substring
   matching, to avoid silently merging unrelated genes.
3. Every combined row is labelled with a "match_type":
     EXACT          - normalised symbols were identical without bla-stripping
     BLA_PREFIX     - symbols only matched after stripping the "bla" prefix
     SINGLE_SOURCE  - gene was called by only one of the two tools
   so you can audit any bla-prefix-assisted matches by eye before trusting them.

Directory conventions (matching PA_comparative.bash)
-----------------------------------------------------
  abricate input  : output of merge_abricate.py   (GENE, DATABASE, ACCESSION,
                     PRODUCT, RESISTANCE, sample1, sample2, ...)
  amrfinder input : output of merge_amrfinder.py  (Element symbol, Element
                     name, Type, Subtype, Class, Subclass, sample1, sample2, ...)

Usage examples
--------------
# Minimal: combine a CARD merge and an AMRFinder merge
python3 merge_abricate_amrfinder.py \\
    --abricate  merged_card.tsv \\
    --amrfinder merged_amr.tsv \\
    --output    merged_card_amrfinder.tsv

# Disable the "bla" prefix harmonisation (exact normalised-name matches only)
python3 merge_abricate_amrfinder.py \\
    --abricate  merged_card.tsv \\
    --amrfinder merged_amr.tsv \\
    --output    merged_card_amrfinder.tsv \\
    --keep-bla-prefix
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

ABRICATE_DESC_COLS = ["GENE", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"]
AMRFINDER_DESC_COLS = ["Element symbol", "Element name", "Type", "Subtype", "Class", "Subclass"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalise(text: str) -> str:
    """Lowercase and strip hyphens, underscores, spaces for fuzzy matching."""
    return re.sub(r"[-_ ]", "", str(text).lower())


def gene_key(symbol: str, strip_bla_prefix: bool) -> str:
    """
    Build the key used to match genes across Abricate and AMRFinder.

    Applies normalise() first, then optionally strips a leading "bla"
    (the standard beta-lactamase gene-name prefix) so e.g. "blaNDM-1"
    (Abricate) and "NDM-1" (AMRFinder) resolve to the same key.
    """
    key = normalise(symbol)
    if strip_bla_prefix and key.startswith("bla"):
        key = key[3:]
    return key


def load_table(path: Path, sep: str, required_cols: list[str], label: str) -> pd.DataFrame:
    """Read a merged wide table and verify it has the expected description columns."""
    if not path.exists():
        sys.exit(f"ERROR: {label} file not found: {path}")

    df = pd.read_csv(path, sep=sep, dtype=str)

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        sys.exit(
            f"ERROR: {label} file ({path}) is missing expected column(s) {missing}. "
            f"Is this really the output of the matching merge script?"
        )
    return df


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def prepare_source(
    df: pd.DataFrame,
    desc_cols: list[str],
    gene_col: str,
    strip_bla_prefix: bool,
) -> tuple[dict[str, dict], dict[str, str], list[str]]:
    """
    Index one source table by its harmonised gene key.

    Returns:
      by_key       - {gene_key: row_dict} (sample counts coerced to int; if two
                      original symbols collide onto the same key, their sample
                      counts are summed and the first-seen description wins)
      raw_norm_key - {gene_key: normalise(original_symbol)} (pre-bla-strip form,
                      used later to label a match as EXACT vs BLA_PREFIX)
      sample_cols  - list of sample columns present in this table
    """
    sample_cols = [c for c in df.columns if c not in desc_cols]

    by_key: dict[str, dict] = {}
    raw_norm_key: dict[str, str] = {}

    for _, row in df.iterrows():
        symbol = row[gene_col]
        if pd.isna(symbol) or str(symbol).strip() == "":
            continue
        key = gene_key(symbol, strip_bla_prefix)
        raw_norm_key.setdefault(key, normalise(symbol))

        if key not in by_key:
            rec = {c: row[c] for c in desc_cols}
            for s in sample_cols:
                rec[s] = int(float(row[s])) if pd.notna(row[s]) and str(row[s]) != "" else 0
            by_key[key] = rec
        else:
            # Rare: two differently-spelled symbols normalised to the same key
            # within one source table. Sum the counts, keep the first symbol.
            for s in sample_cols:
                by_key[key][s] += int(float(row[s])) if pd.notna(row[s]) and str(row[s]) != "" else 0

    return by_key, raw_norm_key, sample_cols


def build_combined_table(
    abricate_by_key: dict[str, dict],
    abricate_raw_norm: dict[str, str],
    abricate_samples: list[str],
    amrfinder_by_key: dict[str, dict],
    amrfinder_raw_norm: dict[str, str],
    amrfinder_samples: list[str],
) -> pd.DataFrame:
    """Combine the two indexed sources into one gene x sample table."""

    # Union of samples, Abricate's ordering first, then any AMRFinder-only samples.
    union_samples = list(abricate_samples)
    for s in amrfinder_samples:
        if s not in union_samples:
            union_samples.append(s)

    only_in_abricate = set(abricate_samples) - set(amrfinder_samples)
    only_in_amrfinder = set(amrfinder_samples) - set(abricate_samples)
    if only_in_abricate:
        print(
            f"WARNING: {len(only_in_abricate)} sample(s) present in the Abricate table "
            f"but not the AMRFinder table (AMRFinder counts will read 0): "
            f"{sorted(only_in_abricate)}",
            file=sys.stderr,
        )
    if only_in_amrfinder:
        print(
            f"WARNING: {len(only_in_amrfinder)} sample(s) present in the AMRFinder table "
            f"but not the Abricate table (Abricate counts will read 0): "
            f"{sorted(only_in_amrfinder)}",
            file=sys.stderr,
        )

    # Preserve first-seen key order: Abricate keys first, then AMRFinder-only keys.
    key_order = list(abricate_by_key.keys())
    for k in amrfinder_by_key:
        if k not in abricate_by_key:
            key_order.append(k)

    records = []
    for key in key_order:
        a_row = abricate_by_key.get(key)
        f_row = amrfinder_by_key.get(key)

        source_parts = []
        if a_row is not None:
            source_parts.append("ABRICATE_CARD")
        if f_row is not None:
            source_parts.append("AMRFINDER")
        source = ";".join(source_parts)

        if a_row is not None and f_row is not None:
            match_type = "EXACT" if abricate_raw_norm[key] == amrfinder_raw_norm[key] else "BLA_PREFIX"
        else:
            match_type = "SINGLE_SOURCE"

        rec = {
            "gene_key": key,
            "abricate_gene": a_row["GENE"] if a_row else "",
            "amrfinder_gene": f_row["Element symbol"] if f_row else "",
            "source": source,
            "match_type": match_type,
            "abricate_database": a_row["DATABASE"] if a_row else "",
            "abricate_accession": a_row["ACCESSION"] if a_row else "",
            "abricate_product": a_row["PRODUCT"] if a_row else "",
            "abricate_resistance": a_row["RESISTANCE"] if a_row else "",
            "amrfinder_name": f_row["Element name"] if f_row else "",
            "amrfinder_type": f_row["Type"] if f_row else "",
            "amrfinder_subtype": f_row["Subtype"] if f_row else "",
            "amrfinder_class": f_row["Class"] if f_row else "",
            "amrfinder_subclass": f_row["Subclass"] if f_row else "",
        }
        for sample in union_samples:
            rec[f"{sample}_abricate"] = a_row[sample] if (a_row and sample in a_row) else 0
            rec[f"{sample}_amrfinder"] = f_row[sample] if (f_row and sample in f_row) else 0

        records.append(rec)

    desc_out_cols = [
        "gene_key", "abricate_gene", "amrfinder_gene", "source", "match_type",
        "abricate_database", "abricate_accession", "abricate_product", "abricate_resistance",
        "amrfinder_name", "amrfinder_type", "amrfinder_subtype", "amrfinder_class", "amrfinder_subclass",
    ]
    sample_out_cols = [f"{s}_{tool}" for s in union_samples for tool in ("abricate", "amrfinder")]

    combined = pd.DataFrame.from_records(records, columns=desc_out_cols + sample_out_cols)
    combined = combined.sort_values("gene_key", kind="stable").reset_index(drop=True)
    return combined


def write_summary(combined: pd.DataFrame, outpath: Path) -> None:
    both = combined[combined["source"] == "ABRICATE_CARD;AMRFINDER"]
    only_abricate = combined[combined["source"] == "ABRICATE_CARD"]
    only_amrfinder = combined[combined["source"] == "AMRFINDER"]
    bla_assisted = combined[combined["match_type"] == "BLA_PREFIX"]

    lines = [
        "=============================================================",
        " Abricate (CARD) vs AMRFinderPlus Reconciliation Summary",
        "=============================================================",
        f"  Total unique genes (by harmonised key) : {len(combined)}",
        f"  Called by BOTH tools                   : {len(both)}",
        f"  Called by ABRICATE_CARD only            : {len(only_abricate)}",
        f"  Called by AMRFINDER only                : {len(only_amrfinder)}",
        f"  Matches only after bla-prefix stripping : {len(bla_assisted)}",
    ]
    if len(bla_assisted) > 0:
        lines.append("")
        lines.append("  bla-prefix-assisted matches (review before trusting):")
        for _, r in bla_assisted.iterrows():
            lines.append(f"    {r['abricate_gene']}  <->  {r['amrfinder_gene']}")
    lines.append("=============================================================")

    print("\n".join(lines))
    summary_path = outpath.with_name(outpath.stem + ".summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  Summary file : {summary_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile a merged Abricate table and a merged AMRFinderPlus table "
            "into one gene x sample table tagged with which tool(s) called each gene."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--abricate", required=True, type=Path, metavar="TSV",
        help="Merged Abricate table, i.e. the --output of merge_abricate.py "
             "(must contain GENE, DATABASE, ACCESSION, PRODUCT, RESISTANCE + sample columns).",
    )
    parser.add_argument(
        "--amrfinder", required=True, type=Path, metavar="TSV",
        help="Merged AMRFinderPlus table, i.e. the --output of merge_amrfinder.py "
             "(must contain Element symbol, Element name, Type, Subtype, Class, Subclass "
             "+ sample columns).",
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="Output combined table path. A companion <output>.summary.txt is also written.",
    )
    parser.add_argument(
        "--abricate-sep", default="\t",
        help="Field separator used in the --abricate input file. Default: tab.",
    )
    parser.add_argument(
        "--amrfinder-sep", default="\t",
        help="Field separator used in the --amrfinder input file. Default: tab.",
    )
    parser.add_argument(
        "--sep", default="\t",
        help="Field separator for the OUTPUT table. Default: tab. Use ',' for CSV.",
    )
    parser.add_argument(
        "--keep-bla-prefix", action="store_true",
        help="Disable the 'bla' prefix harmonisation rule; genes are matched purely on "
             "normalise()'d symbols (lowercased, hyphens/underscores/spaces stripped). "
             "Use this if the bla-prefix rule is producing incorrect matches for your data.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    strip_bla_prefix = not args.keep_bla_prefix

    print(f"\nAbricate table  : {args.abricate}")
    print(f"AMRFinder table : {args.amrfinder}")
    print(f"bla-prefix rule : {'enabled' if strip_bla_prefix else 'disabled (exact match only)'}")
    print(f"Output          : {args.output}\n")

    abricate_df = load_table(args.abricate, args.abricate_sep, ABRICATE_DESC_COLS, "Abricate")
    amrfinder_df = load_table(args.amrfinder, args.amrfinder_sep, AMRFINDER_DESC_COLS, "AMRFinder")

    abricate_by_key, abricate_raw_norm, abricate_samples = prepare_source(
        abricate_df, ABRICATE_DESC_COLS, "GENE", strip_bla_prefix
    )
    amrfinder_by_key, amrfinder_raw_norm, amrfinder_samples = prepare_source(
        amrfinder_df, AMRFINDER_DESC_COLS, "Element symbol", strip_bla_prefix
    )

    combined = build_combined_table(
        abricate_by_key, abricate_raw_norm, abricate_samples,
        amrfinder_by_key, amrfinder_raw_norm, amrfinder_samples,
    )

    if combined.empty:
        sys.exit("ERROR: no genes found in either input table - check the input files.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep=args.sep, index=False)
    print(f"Combined table: {len(combined)} unique genes -> {args.output}\n")

    write_summary(combined, args.output)

    print("\nDone.")


if __name__ == "__main__":
    main()

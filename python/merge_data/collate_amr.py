#!/usr/bin/env python3
"""
collate_amr.py
--------------
Collate per-sample ABRicate and AMRFinderPlus outputs into one reconciled
gene x sample AMR matrix, with optional per-tool copy-count matrices.

BACKGROUND
----------
ABRicate (here run against CARD) and AMRFinderPlus use different reference
databases and different gene-naming conventions for the same resistance gene
(e.g. "blaNDM-1" vs "NDM-1"). Comparing them side by side is the quickest way
to see which resistance calls are supported by both tools (high confidence)
and which come from one tool only. This script replaces the former three-step
chain merge_abricate.py -> merge_amrfinder.py -> merge_abricate_amrfinder.py:
it reads the raw per-sample tool outputs directly, so no intermediate merged
tables are needed (they can still be written with --per-tool-prefix).

GENE-NAME MATCHING STRATEGY
---------------------------
1. normalise(): lowercase, strip hyphens/underscores/spaces.
   "blaNDM-1" -> "blandm1".
2. Unless --keep-bla-prefix is set, a leading "bla" is then stripped so that
   "blandm1" (ABRicate) and "ndm1" (AMRFinderPlus) resolve to the key "ndm1".
   This is a single, narrow, well-understood rule (the beta-lactamase prefix
   convention), not open-ended fuzzy matching.
3. Every reconciled row carries a "match_type":
     EXACT          - normalised symbols identical without bla-stripping
     BLA_PREFIX     - symbols matched only after stripping "bla" (audit these)
     SINGLE_SOURCE  - gene called by only one of the two tools

CELL VALUES
-----------
Reconciled matrix (default): one column per sample, cell = "X,Y" where
  X = ABRicate copy count, Y = AMRFinderPlus copy count.
  "1,1" concordant | "1,0" ABRicate only | "0,1" AMRFinderPlus only | "0,0" absent
  "NA" replaces a count when that tool's file is missing for that sample, so a
  missing file is never silently reported as "gene absent".
With --split-columns, two columns per sample (<sample>_abricate, <sample>_amrfinder).
If only one manifest is supplied, --output receives that tool's plain
copy-count matrix instead of the reconciled matrix.

SAMPLE NAMES & COHORT
---------------------
Manifests may be "sample<TAB>path" (preferred, written by the orchestrator) or
a bare path per line (sample inferred from the file name by stripping known
tool suffixes such as .abricate.card.tsv / .amrfinder.tsv). With --samples,
only the listed isolates are reported, in that order; any other isolate found
in the manifests is dropped and listed in the summary.

Directory conventions (matching Ecoli_merge_script_TN.bash)
-----------------------------------------------------------
  abricate manifest  : {PATHS_DIR}/paths_abricate_card_TN.txt
  amrfinder manifest : {PATHS_DIR}/paths_amrfinder_TN.txt
  cohort sample list : {PATHS_DIR}/samples_TN.txt
  reconciled matrix  : {METADATA_OUT_DIR}/amr_reconciled_TN.tsv (+ .summary.txt)
  per-tool matrices  : {INDIVIDUAL_OUT_DIR}/amr_TN_abricate_matrix.tsv
                       {INDIVIDUAL_OUT_DIR}/amr_TN_amrfinder_matrix.tsv

Usage examples
--------------
# Example 1: reconciled matrix only
python collate_amr.py \\
    --output             merge_metadata/amr_reconciled_TN.tsv \\
    --abricate-manifest  paths/paths_abricate_card_TN.txt \\
    --amrfinder-manifest paths/paths_amrfinder_TN.txt \\
    --samples            paths/samples_TN.txt

# Example 2: also keep per-tool matrices, AMR-type AMRFinder hits only
python collate_amr.py \\
    --output                merge_metadata/amr_reconciled_TN.tsv \\
    --abricate-manifest     paths/paths_abricate_card_TN.txt \\
    --amrfinder-manifest    paths/paths_amrfinder_TN.txt \\
    --samples               paths/samples_TN.txt \\
    --per-tool-prefix       merge_individual_outputs/amr_TN \\
    --amrfinder-type-filter AMR
"""

import argparse
import csv
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# File-name tokens stripped (right to left) when inferring a sample name.
# Identical list in collate_amr.py, collate_genome_overview.py and
# collate_annotations.py so every collator derives the same sample ID.
SAMPLE_SUFFIX_TOKENS = {
    "tsv", "csv", "txt", "gff", "gff3", "gbk", "fasta", "fna", "fa", "faa",
    "annotations", "emapper", "eggnog", "prokka", "bakta",
    "abricate", "card", "vfdb", "ncbi", "resfinder", "plasmidfinder", "amrfinder",
    "quast", "raw", "filtered", "contigs", "scaffolds",
    "isescan", "integron", "integrons", "summary", "phispy",
    "ectyper", "mlst", "kleborate", "checkm",
}

# Generic per-tool file names that carry no sample ID (use the parent directory).
GENERIC_BASENAMES = {"report", "output", "prophage", "prophage_coordinates", "summary", "results"}

ABRICATE_DESC = ["DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"]
AMRFINDER_DESC = ["Element name", "Type", "Subtype", "Class", "Subclass"]


def normalise(text: str) -> str:
    """Lowercase and strip hyphens, underscores and spaces for name matching."""
    return re.sub(r"[-_ ]", "", str(text).lower())


def gene_key(symbol: str, strip_bla_prefix: bool) -> str:
    """
    Build the cross-tool matching key for a gene symbol.

    Args:
        symbol: Raw gene symbol from ABRicate or AMRFinderPlus.
        strip_bla_prefix: If True, drop a leading "bla" after normalisation.

    Returns:
        Harmonised key, e.g. "blaNDM-1" -> "ndm1".
    """
    key = normalise(symbol)
    if strip_bla_prefix and key.startswith("bla"):
        key = key[3:]
    return key


def clean_sample_name(raw: str) -> str:
    """
    Reduce a file name, path or tool label to a bare sample ID.

    Takes the basename, removes an IntegronFinder 'Results_Integron_Finder_'
    prefix, then strips known tool tokens from the right, e.g.
    '132895-LJF30331.prokka.eggnog.emapper.annotations' -> '132895-LJF30331'.

    Args:
        raw: File name, path, or label.

    Returns:
        Cleaned sample ID.
    """
    name = str(raw).strip().replace("\\", "/").rstrip("/").split("/")[-1]
    name = re.sub(r"^Results_Integron_Finder_", "", name, flags=re.IGNORECASE)
    while "." in name:
        stem, _, ext = name.rpartition(".")
        if stem and ext.lower() in SAMPLE_SUFFIX_TOKENS:
            name = stem
        else:
            break
    return name


def sample_from_path(path: Path) -> str:
    """Infer a sample ID from a path, falling back to the parent directory for generic names."""
    name = clean_sample_name(path.name)
    if name.lower() in GENERIC_BASENAMES:
        name = clean_sample_name(path.parent.name)
    return name


def first_present(row: dict, keys: tuple[str, ...]) -> str:
    """Return the first non-empty value among several candidate column names."""
    for k in keys:
        val = row.get(k, "") or ""
        if val.strip():
            return val.strip()
    return ""


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def read_manifest(manifest_path: Path) -> list[tuple[str, Path]]:
    """
    Parse a manifest of 'sample<TAB>path' or bare-path lines.

    Only TAB separates fields, so paths containing spaces are safe.
    Blank lines and lines starting with '#' are ignored.

    Args:
        manifest_path: Manifest text file.

    Returns:
        List of (sample, path) tuples in file order.
    """
    entries: list[tuple[str, Path]] = []
    if not manifest_path.exists():
        print(f"WARNING: manifest not found: {manifest_path}", file=sys.stderr)
        return entries
    with open(manifest_path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t") if p.strip()]
            if len(parts) >= 2:
                entries.append((clean_sample_name(parts[0]), Path(parts[1])))
            else:
                p = Path(parts[0])
                entries.append((sample_from_path(p), p))
    return entries


def read_sample_list(path: Path | None) -> list[str]:
    """Read the cohort sample list (first TAB field of each non-comment line)."""
    if path is None:
        return []
    if not path.exists():
        sys.exit(f"ERROR: --samples file not found: {path}")
    samples: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line and not line.startswith("#"):
                s = clean_sample_name(line.split("\t")[0])
                if s not in samples:
                    samples.append(s)
    return samples


def parse_abricate_tsv(path: Path, database_filter: str | None) -> list[dict] | None:
    """
    Parse one ABRicate per-sample TSV (default 15-column format).

    ABRicate columns (tab-separated):
        #FILE  SEQUENCE  START  END  STRAND  GENE  COVERAGE  COVERAGE_MAP  GAPS
        %COVERAGE  %IDENTITY  DATABASE  ACCESSION  PRODUCT  RESISTANCE

    Args:
        path: Per-sample ABRicate TSV.
        database_filter: Keep only rows whose DATABASE matches (case-insensitive).

    Returns:
        List of {"symbol", DESC...} dicts (one per hit), or None if unreadable.
    """
    if not path.is_file():
        return None
    hits: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = [f.lstrip("#") for f in (reader.fieldnames or [])]
        reader.fieldnames = fields
        if "GENE" not in fields:
            print(f"WARNING: {path} has no GENE column - not an ABRicate TSV, skipping",
                  file=sys.stderr)
            return None
        for row in reader:
            symbol = (row.get("GENE", "") or "").strip()
            if not symbol:
                continue
            if database_filter and (row.get("DATABASE", "") or "").lower() != database_filter.lower():
                continue
            rec = {"symbol": symbol}
            for col in ABRICATE_DESC:
                rec[col] = row.get(col, "") or ""
            hits.append(rec)
    return hits


def parse_amrfinder_tsv(path: Path, type_filter: list[str] | None) -> list[dict] | None:
    """
    Parse one AMRFinderPlus per-sample TSV (v3 'Gene symbol' or v4 'Element symbol').

    AMRFinderPlus columns used (tab-separated):
        Element symbol | Gene symbol, Element name | Sequence name,
        Type, Subtype, Class, Subclass

    Args:
        path: Per-sample AMRFinderPlus TSV.
        type_filter: Keep only rows whose Type is in this list (e.g. ["AMR"]).

    Returns:
        List of {"symbol", DESC...} dicts (one per hit), or None if unreadable.
    """
    if not path.is_file():
        return None
    hits: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        if "Element symbol" not in fields and "Gene symbol" not in fields:
            print(f"WARNING: {path} has no 'Element symbol'/'Gene symbol' column - skipping",
                  file=sys.stderr)
            return None
        for row in reader:
            symbol = first_present(row, ("Element symbol", "Gene symbol"))
            if not symbol:
                continue
            if type_filter and (row.get("Type", "") or "") not in type_filter:
                continue
            hits.append({
                "symbol": symbol,
                "Element name": first_present(row, ("Element name", "Sequence name")),
                "Type": row.get("Type", "") or "",
                "Subtype": row.get("Subtype", "") or "",
                "Class": row.get("Class", "") or "",
                "Subclass": row.get("Subclass", "") or "",
            })
    return hits


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def load_tool(
    entries: list[tuple[str, Path]],
    tool: str,
    database_filter: str | None,
    type_filter: list[str] | None,
) -> tuple[dict[str, dict], dict[str, dict[str, int]]]:
    """
    Load every per-sample file of one tool into description and count tables.

    Args:
        entries: (sample, path) tuples from the manifest.
        tool: "abricate" or "amrfinder".
        database_filter: ABRicate DATABASE filter (ignored for AMRFinderPlus).
        type_filter: AMRFinderPlus Type filter (ignored for ABRicate).

    Returns:
        (desc_by_symbol, counts_by_sample) where counts_by_sample only contains
        samples whose file was read successfully.
    """
    desc_by_symbol: dict[str, dict] = {}
    counts_by_sample: dict[str, dict[str, int]] = {}
    for sample, path in entries:
        if tool == "abricate":
            hits = parse_abricate_tsv(path, database_filter)
        else:
            hits = parse_amrfinder_tsv(path, type_filter)
        if hits is None:
            print(f"  [{sample}] WARNING: {tool} file missing/unreadable: {path} - cells set to NA",
                  file=sys.stderr)
            continue
        if sample in counts_by_sample:
            print(f"  [{sample}] WARNING: duplicate {tool} entry - later file overwrites earlier",
                  file=sys.stderr)
        counts: dict[str, int] = {}
        for h in hits:
            counts[h["symbol"]] = counts.get(h["symbol"], 0) + 1
            desc_by_symbol.setdefault(h["symbol"], h)
        counts_by_sample[sample] = counts
    return desc_by_symbol, counts_by_sample


def resolve_cohort(
    sample_list: list[str],
    manifest_samples: list[str],
) -> tuple[list[str], list[str]]:
    """
    Decide the reported sample columns and the samples dropped by the cohort filter.

    Args:
        sample_list: Samples from --samples (may be empty).
        manifest_samples: Samples seen in the manifests, in first-seen order.

    Returns:
        (reported_samples, dropped_samples).
    """
    if not sample_list:
        return manifest_samples, []
    keep = set(sample_list)
    dropped = [s for s in manifest_samples if s not in keep]
    return sample_list, dropped


def build_tool_matrix(
    desc_by_symbol: dict[str, dict],
    counts_by_sample: dict[str, dict[str, int]],
    samples: list[str],
    desc_cols: list[str],
    symbol_col: str,
) -> list[dict]:
    """
    Build a plain gene x sample copy-count matrix for a single tool.

    Returns:
        One dict per gene symbol (sorted), sample cells are int or "NA".
    """
    rows: list[dict] = []
    for symbol in sorted(desc_by_symbol, key=str.lower):
        desc = desc_by_symbol[symbol]
        rec = {symbol_col: symbol}
        for col in desc_cols:
            rec[col] = desc.get(col, "")
        present = False
        for s in samples:
            if s in counts_by_sample:
                n = counts_by_sample[s].get(symbol, 0)
                rec[s] = n
                present = present or n > 0
            else:
                rec[s] = "NA"
        if present:
            rows.append(rec)
    return rows


def index_by_key(
    desc_by_symbol: dict[str, dict],
    counts_by_sample: dict[str, dict[str, int]],
    strip_bla_prefix: bool,
) -> dict[str, dict]:
    """
    Re-index one tool's results by harmonised gene key.

    If two spellings collapse onto one key within the same tool, their counts
    are summed and the first-seen symbol/description is kept.

    Returns:
        {key: {"symbol", "raw_norm", "desc", "counts": {sample: int}}}
    """
    by_key: dict[str, dict] = {}
    for symbol in desc_by_symbol:
        key = gene_key(symbol, strip_bla_prefix)
        entry = by_key.setdefault(key, {
            "symbol": symbol,
            "raw_norm": normalise(symbol),
            "desc": desc_by_symbol[symbol],
            "counts": {},
        })
        for sample, counts in counts_by_sample.items():
            n = counts.get(symbol, 0)
            if n:
                entry["counts"][sample] = entry["counts"].get(sample, 0) + n
    return by_key


def build_reconciled_matrix(
    abr_keys: dict[str, dict],
    amr_keys: dict[str, dict],
    abr_loaded: set[str],
    amr_loaded: set[str],
    samples: list[str],
    split_columns: bool,
) -> list[dict]:
    """
    Combine ABRicate and AMRFinderPlus key indexes into one reconciled matrix.

    Returns:
        One dict per harmonised gene key that is present in >= 1 reported sample.
    """
    rows: list[dict] = []
    for key in sorted(set(abr_keys) | set(amr_keys)):
        a = abr_keys.get(key)
        f = amr_keys.get(key)
        if a and f:
            source = "ABRICATE_CARD;AMRFINDER"
            match_type = "EXACT" if a["raw_norm"] == f["raw_norm"] else "BLA_PREFIX"
        else:
            source = "ABRICATE_CARD" if a else "AMRFINDER"
            match_type = "SINGLE_SOURCE"

        a_desc = a["desc"] if a else {}
        f_desc = f["desc"] if f else {}
        rec = {
            "gene_key": key,
            "abricate_gene": a["symbol"] if a else "",
            "amrfinder_gene": f["symbol"] if f else "",
            "source": source,
            "match_type": match_type,
            "abricate_database": a_desc.get("DATABASE", ""),
            "abricate_accession": a_desc.get("ACCESSION", ""),
            "abricate_product": a_desc.get("PRODUCT", ""),
            "abricate_resistance": a_desc.get("RESISTANCE", ""),
            "amrfinder_name": f_desc.get("Element name", ""),
            "amrfinder_type": f_desc.get("Type", ""),
            "amrfinder_subtype": f_desc.get("Subtype", ""),
            "amrfinder_class": f_desc.get("Class", ""),
            "amrfinder_subclass": f_desc.get("Subclass", ""),
        }

        n_abr = n_amr = n_both = 0
        for s in samples:
            a_cnt = (a["counts"].get(s, 0) if a else 0) if s in abr_loaded else "NA"
            f_cnt = (f["counts"].get(s, 0) if f else 0) if s in amr_loaded else "NA"
            a_pos = isinstance(a_cnt, int) and a_cnt > 0
            f_pos = isinstance(f_cnt, int) and f_cnt > 0
            n_abr += a_pos
            n_amr += f_pos
            n_both += a_pos and f_pos
            if split_columns:
                rec[f"{s}_abricate"] = a_cnt
                rec[f"{s}_amrfinder"] = f_cnt
            else:
                rec[s] = f"{a_cnt},{f_cnt}"
        rec["n_samples_abricate"] = n_abr
        rec["n_samples_amrfinder"] = n_amr
        rec["n_samples_concordant"] = n_both
        if n_abr or n_amr:
            rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

RECONCILED_DESC_COLUMNS = [
    "gene_key", "abricate_gene", "amrfinder_gene", "source", "match_type",
    "n_samples_abricate", "n_samples_amrfinder", "n_samples_concordant",
    "abricate_database", "abricate_accession", "abricate_product", "abricate_resistance",
    "amrfinder_name", "amrfinder_type", "amrfinder_subtype", "amrfinder_class", "amrfinder_subclass",
]


def write_table(
    rows: list[dict],
    columns: list[str],
    outpath: Path,
    comments: list[str] | None = None,
) -> None:
    """Write rows to a TSV with optional leading '# ' comment lines."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        for c in comments or []:
            fh.write(f"# {c}\n")
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {outpath} ({len(rows)} rows x {len(columns)} columns)")


def write_summary(
    rows: list[dict],
    samples: list[str],
    abr_loaded: set[str],
    amr_loaded: set[str],
    dropped: list[str],
    strip_bla_prefix: bool,
    outpath: Path,
) -> None:
    """Write and echo the reconciliation audit summary."""
    both = [r for r in rows if r["source"] == "ABRICATE_CARD;AMRFINDER"]
    only_a = [r for r in rows if r["source"] == "ABRICATE_CARD"]
    only_f = [r for r in rows if r["source"] == "AMRFINDER"]
    bla = [r for r in rows if r["match_type"] == "BLA_PREFIX"]
    missing_a = [s for s in samples if s not in abr_loaded]
    missing_f = [s for s in samples if s not in amr_loaded]

    lines = [
        "=============================================================",
        " ABRicate (CARD) vs AMRFinderPlus Reconciliation Summary",
        "=============================================================",
        f"  Reported samples                       : {len(samples)}",
        f"  Samples missing ABRicate file (NA)     : {len(missing_a)} {missing_a if missing_a else ''}",
        f"  Samples missing AMRFinderPlus file (NA): {len(missing_f)} {missing_f if missing_f else ''}",
        f"  Samples dropped by --samples filter    : {len(dropped)} {dropped if dropped else ''}",
        f"  bla-prefix rule                        : {'enabled' if strip_bla_prefix else 'disabled'}",
        "-------------------------------------------------------------",
        f"  Unique genes (harmonised key)          : {len(rows)}",
        f"  Called by BOTH tools                   : {len(both)}",
        f"  Called by ABRICATE_CARD only           : {len(only_a)}",
        f"  Called by AMRFINDER only               : {len(only_f)}",
        f"  Matched only after bla-prefix strip    : {len(bla)}",
    ]
    if bla:
        lines.append("  bla-prefix-assisted matches (review before trusting):")
        lines.extend(f"    {r['abricate_gene']}  <->  {r['amrfinder_gene']}" for r in bla)
    lines.append("=============================================================")

    with open(outpath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n  Summary file : {outpath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Collate per-sample ABRicate and AMRFinderPlus outputs into one reconciled "
            "gene x sample AMR matrix, with optional per-tool copy-count matrices."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Required arguments ---
    parser.add_argument(
        "--output", required=True, type=Path, metavar="TSV",
        help=(
            "Destination of the main matrix. With both manifests this is the reconciled "
            "'X,Y' matrix; with one manifest it is that tool's copy-count matrix. "
            "A companion <output>.summary.txt is written next to it "
            "(orchestrator: ${METADATA_OUT_DIR}/amr_reconciled_TN.tsv)."
        ),
    )

    # --- Optional inputs ---
    parser.add_argument(
        "--abricate-manifest", type=Path, default=None, metavar="TXT",
        help=(
            "Manifest of per-sample ABRicate TSVs ('sample<TAB>path' or bare path per line). "
            "Point it at one database's files (e.g. every *.abricate.card.tsv). "
            "Orchestrator: ${PATHS_DIR}/paths_abricate_card_TN.txt."
        ),
    )
    parser.add_argument(
        "--amrfinder-manifest", type=Path, default=None, metavar="TXT",
        help=(
            "Manifest of per-sample AMRFinderPlus TSVs ('sample<TAB>path' or bare path). "
            "Orchestrator: ${PATHS_DIR}/paths_amrfinder_TN.txt."
        ),
    )
    parser.add_argument(
        "--samples", type=Path, default=None, metavar="TXT",
        help=(
            "Optional cohort list (one sample ID per line). Only these isolates are "
            "reported, in this order; others found in the manifests are dropped and "
            "listed in the summary. Orchestrator: ${PATHS_DIR}/samples_TN.txt."
        ),
    )
    parser.add_argument(
        "--per-tool-prefix", type=Path, default=None, metavar="PREFIX",
        help=(
            "If set (and both manifests are given), also write the per-tool copy-count "
            "matrices to <PREFIX>_abricate_matrix.tsv and <PREFIX>_amrfinder_matrix.tsv. "
            "Orchestrator: ${INDIVIDUAL_OUT_DIR}/amr_TN."
        ),
    )

    # --- Filters & layout ---
    parser.add_argument(
        "--abricate-database", default=None, metavar="DB",
        help=(
            "Safety filter: keep only ABRicate rows whose DATABASE column equals this "
            "value (case-insensitive), e.g. 'card'. Default: keep all rows."
        ),
    )
    parser.add_argument(
        "--amrfinder-type-filter", nargs="+", default=None, metavar="TYPE",
        help=(
            "Keep only AMRFinderPlus rows whose Type is one of these values "
            "(e.g. AMR, or AMR STRESS VIRULENCE). Default: keep all types."
        ),
    )
    parser.add_argument(
        "--keep-bla-prefix", action="store_true",
        help=(
            "Disable the 'bla' prefix harmonisation rule; genes then match only on "
            "normalised symbols. Use if bla-assisted matches look wrong."
        ),
    )
    parser.add_argument(
        "--split-columns", action="store_true",
        help=(
            "Write two columns per sample (<sample>_abricate, <sample>_amrfinder) instead "
            "of the default single 'X,Y' column per sample."
        ),
    )

    return parser.parse_args()


def main():
    """Orchestrate loading, reconciliation and writing."""
    args = parse_args()

    if not args.abricate_manifest and not args.amrfinder_manifest:
        sys.exit("ERROR: provide --abricate-manifest and/or --amrfinder-manifest.")

    strip_bla_prefix = not args.keep_bla_prefix
    sample_list = read_sample_list(args.samples)

    print("\n========================================================")
    print(" collate_amr.py - ABRicate + AMRFinderPlus collation")
    print("========================================================")
    print(f" ABRicate manifest  : {args.abricate_manifest or 'None'}")
    print(f" AMRFinder manifest : {args.amrfinder_manifest or 'None'}")
    print(f" Cohort list        : {args.samples or 'None (all manifest samples)'}")
    print(f" bla-prefix rule    : {'enabled' if strip_bla_prefix else 'disabled'}")
    print(f" Layout             : {'2 columns per sample' if args.split_columns else 'X,Y per sample'}")
    print(f" Output             : {args.output}")
    print("========================================================\n")

    abr_entries = read_manifest(args.abricate_manifest) if args.abricate_manifest else []
    amr_entries = read_manifest(args.amrfinder_manifest) if args.amrfinder_manifest else []

    abr_desc, abr_counts = load_tool(abr_entries, "abricate", args.abricate_database, None)
    amr_desc, amr_counts = load_tool(amr_entries, "amrfinder", None, args.amrfinder_type_filter)

    manifest_samples: list[str] = []
    for s, _ in abr_entries + amr_entries:
        if s not in manifest_samples:
            manifest_samples.append(s)
    samples, dropped = resolve_cohort(sample_list, manifest_samples)
    if dropped:
        print(f"Cohort filter dropped {len(dropped)} sample(s): {dropped}", file=sys.stderr)
    if not samples:
        sys.exit("ERROR: no samples to report - check manifests and --samples.")

    # Single-tool mode: write that tool's plain matrix and stop.
    if not (abr_entries and amr_entries):
        if abr_entries:
            rows = build_tool_matrix(abr_desc, abr_counts, samples, ABRICATE_DESC, "GENE")
            cols = ["GENE"] + ABRICATE_DESC + samples
            label = "ABRicate"
        else:
            rows = build_tool_matrix(amr_desc, amr_counts, samples, AMRFINDER_DESC, "Element symbol")
            cols = ["Element symbol"] + AMRFINDER_DESC + samples
            label = "AMRFinderPlus"
        write_table(rows, cols, args.output, [
            f"TOOL: collate_amr.py ({label} copy-count matrix; single-tool mode)",
            "FORMAT: cell = number of hits of the gene in the sample; NA = tool output missing for that sample.",
        ])
        print("\nDone.")
        return

    if args.per_tool_prefix:
        write_table(
            build_tool_matrix(abr_desc, abr_counts, samples, ABRICATE_DESC, "GENE"),
            ["GENE"] + ABRICATE_DESC + samples,
            Path(f"{args.per_tool_prefix}_abricate_matrix.tsv"),
            ["TOOL: collate_amr.py (ABRicate copy-count matrix)",
             "FORMAT: cell = number of hits; NA = ABRicate output missing for that sample."],
        )
        write_table(
            build_tool_matrix(amr_desc, amr_counts, samples, AMRFINDER_DESC, "Element symbol"),
            ["Element symbol"] + AMRFINDER_DESC + samples,
            Path(f"{args.per_tool_prefix}_amrfinder_matrix.tsv"),
            ["TOOL: collate_amr.py (AMRFinderPlus copy-count matrix)",
             "FORMAT: cell = number of hits; NA = AMRFinderPlus output missing for that sample."],
        )

    rows = build_reconciled_matrix(
        abr_keys=index_by_key(abr_desc, abr_counts, strip_bla_prefix),
        amr_keys=index_by_key(amr_desc, amr_counts, strip_bla_prefix),
        abr_loaded=set(abr_counts),
        amr_loaded=set(amr_counts),
        samples=samples,
        split_columns=args.split_columns,
    )
    if not rows:
        sys.exit("ERROR: no genes found in either tool - check the input files.")

    if args.split_columns:
        sample_cols = [f"{s}_{t}" for s in samples for t in ("abricate", "amrfinder")]
        fmt = "FORMAT: 2 columns per sample (<sample>_abricate, <sample>_amrfinder) holding copy counts; NA = that tool's output missing."
    else:
        sample_cols = list(samples)
        fmt = "FORMAT: 1 column per sample, cell = 'X,Y' (X = ABRicate CARD copies, Y = AMRFinderPlus copies); NA = that tool's output missing."
    write_table(rows, RECONCILED_DESC_COLUMNS + sample_cols, args.output, [
        "TOOL: collate_amr.py (ABRicate CARD vs AMRFinderPlus reconciliation matrix)",
        fmt,
        "INTERPRETATION: '1,1' concordant (high confidence); '1,0' ABRicate only; '0,1' AMRFinderPlus only; '0,0' absent.",
        f"MATCHING: normalised symbols{' + leading bla stripped' if strip_bla_prefix else ''}; match_type EXACT / BLA_PREFIX / SINGLE_SOURCE.",
        "COUNTS: n_samples_* = number of reported samples with >= 1 copy by that tool / by both tools.",
    ])

    write_summary(
        rows, samples, set(abr_counts), set(amr_counts), dropped, strip_bla_prefix,
        args.output.with_name(args.output.stem + ".summary.txt"),
    )
    print("\nDone.")


if __name__ == "__main__":
    main()

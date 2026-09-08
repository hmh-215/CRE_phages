#!/usr/bin/env python3
"""
merge_eggnog.py
----------------
Concatenate eggNOG-mapper ".emapper.annotations" tables from multiple
samples into one long table, with a "Sample" column added so every row
is traceable back to its source sample. Optionally also filter that table
down and/or produce compact per-sample summaries, since 50 samples x a
few thousand genes each adds up to a very large concatenated file fast.

INPUT
-----
A plain text "manifest" file listing one sample per line, in either form:

    /path/to/WS2762512A05.eggnog/WS2762512A05.eggnog.emapper.annotations
    WS2762512A05<TAB>/path/to/WS2762512A05.eggnog

The path can point to either:
  - the .emapper.annotations FILE directly, or
  - the eggNOG-mapper output DIRECTORY (the script globs for
    "*.emapper.annotations" inside it and uses the single match found;
    if more than one file matches, the first alphabetically is used and
    a warning is printed).

- If only a path is given, the sample name is inferred from the filename
  or directory name (common suffixes like ".eggnog.emapper.annotations"
  are stripped automatically).
- Blank lines and lines starting with '#' in the MANIFEST are ignored
  (this is separate from the '##' comment lines inside each annotations
  file itself, which are always stripped during parsing).
- Manifest fields can be separated by a tab or by whitespace.

eggNOG-mapper's annotation files contain '##' comment lines (run info,
timestamps, totals) around a single header line that starts with a lone
'#' (e.g. "#query\\tseed_ortholog\\t..."). Both are handled automatically:
'##' lines are dropped, and the leading '#' is stripped from the header
so the first column comes through as a clean "query" column.

OUTPUTS (pick any combination - at least one is required)
-----------------------------------------------------------
--full-output PATH
    The concatenated per-gene table (optionally filtered/trimmed - see
    below). Skip this flag entirely if you only want a summary and don't
    need the raw rows.

--summary cog kegg completeness   (space-separated, choose one or more)
--summary-output-prefix PREFIX
    Writes PREFIX_cog.tsv / PREFIX_kegg.tsv / PREFIX_completeness.tsv as
    requested. Summaries are always computed from the FULL, unfiltered
    data (row filters below only affect --full-output), so percentages
    stay meaningful regardless of how you trim the raw table.

    cog          : per-sample counts of each COG functional category
                   letter (multi-letter COG_category strings like "KT"
                   are split and each letter counted once per gene).
    kegg         : per-sample counts of genes assigned to each KEGG
                   pathway (from the KEGG_Pathway column; a gene listed
                   under several pathways is counted once per pathway).
    completeness : one row per sample with total gene count and the
                   number/percent of genes that have a hit, a COG
                   category, a GO term, a KEGG ortholog, and a PFAM hit.

FULL-TABLE ROW FILTERS (only apply to --full-output)
-----------------------------------------------------
--exclude-unannotated       drop rows with no eggNOG hit at all (Description == "-")
--min-score FLOAT           keep only rows with 'score' >= this value
--max-evalue FLOAT          keep only rows with 'evalue' <= this value
--contains TEXT             keep only rows where TEXT appears (case-insensitive,
                             regex allowed) in query, Description, or Preferred_name
--keep-columns COL [COL...] trim the full output to just these columns
                             ("Sample" is always kept even if you omit it)

USAGE
-----
    # just the summaries, skip the giant raw table entirely
    python3 merge_eggnog.py --manifest samples.txt \\
        --summary cog completeness --summary-output-prefix results/eggnog

    # full table, but trimmed to annotated hits only and fewer columns
    python3 merge_eggnog.py --manifest samples.txt \\
        --full-output merged_eggnog.filtered.tsv \\
        --exclude-unannotated \\
        --keep-columns query Description Preferred_name COG_category GOs KEGG_Pathway

    # search for a gene/function of interest across all 50 samples
    python3 merge_eggnog.py --manifest samples.txt \\
        --full-output quorum_sensing_hits.tsv --contains "quorum sensing"

    # everything at once
    python3 merge_eggnog.py --manifest samples.txt \\
        --full-output merged_eggnog.tsv \\
        --summary cog kegg completeness --summary-output-prefix results/eggnog
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

NAME_SUFFIXES = [
    ".eggnog.emapper.annotations",
    "_eggnog.emapper.annotations",
    ".emapper.annotations",
    ".eggnog",
    "_eggnog",
]


def infer_sample_name(path: Path) -> str:
    name = path.name
    for suf in NAME_SUFFIXES:
        if name.endswith(suf):
            return name[: -len(suf)]
    return path.stem if path.suffix else path.name


def parse_manifest(manifest_path: Path):
    """Return list of (sample_name, path) tuples. path may be a dir or an annotations file."""
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


def resolve_annotations_file(path: Path) -> Path | None:
    """Given a file or directory path, return the .emapper.annotations file path, or None."""
    if path.is_file():
        return path
    if path.is_dir():
        matches = sorted(path.glob("*.emapper.annotations"))
        if not matches:
            return None
        if len(matches) > 1:
            print(f"WARNING: multiple *.emapper.annotations files found in {path}, "
                  f"using {matches[0].name} and ignoring: "
                  f"{[m.name for m in matches[1:]]}", file=sys.stderr)
        return matches[0]
    return None


def load_annotations(sample: str, ann_path: Path) -> pd.DataFrame:
    """Parse one .emapper.annotations file, stripping '##' comments and cleaning the header."""
    with open(ann_path) as fh:
        raw_lines = [ln.rstrip("\n") for ln in fh]

    # drop eggNOG-mapper's run-info / timestamp / totals comment lines ('##...')
    # but keep the single-'#' header line and all data lines
    kept = [ln for ln in raw_lines if not ln.startswith("##")]
    kept = [ln for ln in kept if ln.strip() != ""]

    if not kept:
        print(f"WARNING: no data found in {ann_path} for sample '{sample}' - skipping",
              file=sys.stderr)
        return pd.DataFrame()

    header_line = kept[0]
    if not header_line.startswith("#"):
        print(f"WARNING: unexpected format in {ann_path} (no leading '#' on header line) "
              f"- skipping sample '{sample}'", file=sys.stderr)
        return pd.DataFrame()

    columns = header_line.lstrip("#").split("\t")
    data_lines = kept[1:]

    rows = [ln.split("\t") for ln in data_lines]
    # guard against any stray short/long rows rather than crashing the whole merge
    rows = [r for r in rows if len(r) == len(columns)]
    if len(rows) < len(data_lines):
        print(f"WARNING: {len(data_lines) - len(rows)} malformed row(s) skipped in {ann_path}",
              file=sys.stderr)

    df = pd.DataFrame(rows, columns=columns)
    df.insert(0, "Sample", sample)
    return df


def build_merged_table(entries) -> pd.DataFrame:
    frames = []
    for sample, path in entries:
        ann_path = resolve_annotations_file(path)
        if ann_path is None:
            print(f"WARNING: could not find a .emapper.annotations file for sample "
                  f"'{sample}' at {path} - skipping", file=sys.stderr)
            continue
        df = load_annotations(sample, ann_path)
        if not df.empty:
            frames.append(df)

    if not frames:
        sys.exit("ERROR: no annotation data could be loaded from any manifest entry.")

    # outer-join columns in case of any version drift between runs, preserving
    # first-seen column order
    merged = pd.concat(frames, ignore_index=True, sort=False)
    return merged


def apply_full_output_filters(df: pd.DataFrame, args) -> pd.DataFrame:
    """Apply the optional row filters + column trim, used only for --full-output."""
    out = df

    if args.exclude_unannotated:
        if "Description" in out.columns:
            before = len(out)
            out = out[out["Description"] != "-"]
            print(f"--exclude-unannotated: dropped {before - len(out)} unannotated row(s)")
        else:
            print("WARNING: --exclude-unannotated requested but no 'Description' column found "
                  "- skipping this filter", file=sys.stderr)

    if args.min_score is not None:
        if "score" in out.columns:
            before = len(out)
            numeric_score = pd.to_numeric(out["score"], errors="coerce")
            out = out[numeric_score >= args.min_score]
            print(f"--min-score {args.min_score}: dropped {before - len(out)} row(s)")
        else:
            print("WARNING: --min-score requested but no 'score' column found - skipping",
                  file=sys.stderr)

    if args.max_evalue is not None:
        if "evalue" in out.columns:
            before = len(out)
            numeric_eval = pd.to_numeric(out["evalue"], errors="coerce")
            out = out[numeric_eval <= args.max_evalue]
            print(f"--max-evalue {args.max_evalue}: dropped {before - len(out)} row(s)")
        else:
            print("WARNING: --max-evalue requested but no 'evalue' column found - skipping",
                  file=sys.stderr)

    if args.contains:
        search_cols = [c for c in ("query", "Description", "Preferred_name") if c in out.columns]
        if search_cols:
            mask = pd.Series(False, index=out.index)
            for c in search_cols:
                mask = mask | out[c].str.contains(args.contains, case=False, regex=True, na=False)
            before = len(out)
            out = out[mask]
            print(f"--contains {args.contains!r}: kept {len(out)} of {before} row(s) "
                  f"(searched columns: {search_cols})")
        else:
            print("WARNING: --contains requested but none of query/Description/Preferred_name "
                  "found - skipping this filter", file=sys.stderr)

    if args.keep_columns:
        missing = [c for c in args.keep_columns if c not in out.columns]
        if missing:
            print(f"WARNING: --keep-columns not found and will be skipped: {missing}",
                  file=sys.stderr)
        keep = ["Sample"] + [c for c in args.keep_columns if c in out.columns and c != "Sample"]
        out = out[keep]

    return out


def summarize_cog(df: pd.DataFrame) -> pd.DataFrame:
    """Per-sample counts of each COG functional category letter."""
    if "COG_category" not in df.columns:
        sys.exit("ERROR: cannot build COG summary - no 'COG_category' column in the data.")

    counts: dict[str, dict[str, int]] = {}
    for sample, cog_str in zip(df["Sample"], df["COG_category"]):
        sample_counts = counts.setdefault(sample, {})
        cog_str = (cog_str or "").strip()
        if cog_str in ("", "-"):
            sample_counts["Unclassified"] = sample_counts.get("Unclassified", 0) + 1
            continue
        for letter in cog_str:
            sample_counts[letter] = sample_counts.get(letter, 0) + 1

    result = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    # sort columns alphabetically but keep "Unclassified" last
    cols = sorted(c for c in result.columns if c != "Unclassified")
    if "Unclassified" in result.columns:
        cols = cols + ["Unclassified"]
    result = result[cols]
    result.index.name = "Sample"
    return result.reset_index()


def summarize_kegg(df: pd.DataFrame) -> pd.DataFrame:
    """Per-sample counts of genes assigned to each KEGG pathway."""
    if "KEGG_Pathway" not in df.columns:
        sys.exit("ERROR: cannot build KEGG summary - no 'KEGG_Pathway' column in the data.")

    counts: dict[str, dict[str, int]] = {}
    for sample, pathway_str in zip(df["Sample"], df["KEGG_Pathway"]):
        sample_counts = counts.setdefault(sample, {})
        pathway_str = (pathway_str or "").strip()
        if pathway_str in ("", "-"):
            sample_counts["Unassigned"] = sample_counts.get("Unassigned", 0) + 1
            continue
        # KEGG_Pathway is a comma-separated list of pathway/module IDs
        for pathway in pathway_str.split(","):
            pathway = pathway.strip()
            if pathway:
                sample_counts[pathway] = sample_counts.get(pathway, 0) + 1

    result = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    cols = sorted(c for c in result.columns if c != "Unassigned")
    if "Unassigned" in result.columns:
        cols = cols + ["Unassigned"]
    result = result[cols]
    result.index.name = "Sample"
    return result.reset_index()


# columns checked for "does this gene have an annotation here" in the completeness summary
COMPLETENESS_FIELDS = [
    ("seed_ortholog", "Any_hit"),
    ("COG_category", "COG_category"),
    ("GOs", "GO_terms"),
    ("KEGG_ko", "KEGG_ortholog"),
    ("PFAMs", "PFAM"),
]


def summarize_completeness(df: pd.DataFrame) -> pd.DataFrame:
    """One row per sample: total genes + count/percent annotated for each key field."""
    available_fields = [(col, label) for col, label in COMPLETENESS_FIELDS if col in df.columns]
    if not available_fields:
        sys.exit("ERROR: cannot build completeness summary - none of the expected columns "
                  f"({[c for c, _ in COMPLETENESS_FIELDS]}) were found in the data.")

    rows = []
    for sample, group in df.groupby("Sample"):
        total = len(group)
        rec = {"Sample": sample, "Total_genes": total}
        for col, label in available_fields:
            n_annotated = int((group[col] != "-").sum())
            pct = round(100 * n_annotated / total, 2) if total else 0.0
            rec[f"{label}_n"] = n_annotated
            rec[f"{label}_pct"] = pct
        rows.append(rec)

    return pd.DataFrame(rows).sort_values("Sample").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(
        description="Concatenate eggNOG-mapper .emapper.annotations files from multiple "
                    "samples into one long table, with optional filtering and per-sample "
                    "summaries (COG categories, KEGG pathways, annotation completeness).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--manifest", required=True, type=Path,
                     help="Text file listing sample annotation files/dirs (one per line, "
                          "optionally 'sample<TAB>path').")
    ap.add_argument("--sep", default="\t",
                     help="Field separator for output tables (default: tab). Use ',' for CSV.")

    # full concatenated table (optional)
    ap.add_argument("--full-output", type=Path, default=None,
                     help="Write the concatenated (optionally filtered/trimmed) per-gene table here. "
                          "Omit this flag if you only want summaries.")
    ap.add_argument("--exclude-unannotated", action="store_true",
                     help="Full-table only: drop rows with no eggNOG hit (Description == '-').")
    ap.add_argument("--min-score", type=float, default=None,
                     help="Full-table only: keep rows with 'score' >= this value.")
    ap.add_argument("--max-evalue", type=float, default=None,
                     help="Full-table only: keep rows with 'evalue' <= this value.")
    ap.add_argument("--contains", type=str, default=None,
                     help="Full-table only: keep rows where this text/regex appears in "
                          "query, Description, or Preferred_name (case-insensitive).")
    ap.add_argument("--keep-columns", nargs="+", default=None,
                     help="Full-table only: trim output to just these columns "
                          "('Sample' is always kept).")

    # summaries (optional)
    ap.add_argument("--summary", nargs="+", choices=["cog", "kegg", "completeness"], default=None,
                     help="One or more summary types to generate (always computed from the "
                          "full, unfiltered data).")
    ap.add_argument("--summary-output-prefix", type=Path, default=None,
                     help="Prefix for summary output files, e.g. 'results/eggnog' produces "
                          "results/eggnog_cog.tsv, results/eggnog_completeness.tsv, etc. "
                          "Required if --summary is given.")

    args = ap.parse_args()

    if not args.full_output and not args.summary:
        ap.error("nothing to do - specify --full-output and/or --summary.")
    if args.summary and not args.summary_output_prefix:
        ap.error("--summary requires --summary-output-prefix.")

    if not args.manifest.exists():
        sys.exit(f"ERROR: manifest file not found: {args.manifest}")

    entries = parse_manifest(args.manifest)
    if not entries:
        sys.exit(f"ERROR: no valid sample entries parsed from {args.manifest}")

    print(f"Parsed {len(entries)} sample entries from manifest.")

    merged = build_merged_table(entries)
    n_rows = len(merged)
    n_samples = merged["Sample"].nunique()
    print(f"Loaded {n_rows} annotation rows from {n_samples} samples.")

    if args.full_output:
        filtered = apply_full_output_filters(merged, args)
        args.full_output.parent.mkdir(parents=True, exist_ok=True)
        filtered.to_csv(args.full_output, sep=args.sep, index=False)
        print(f"Full table: {len(filtered)} rows -> {args.full_output}")

    if args.summary:
        args.summary_output_prefix.parent.mkdir(parents=True, exist_ok=True)
        ext = "csv" if args.sep == "," else "tsv"

        if "cog" in args.summary:
            cog_df = summarize_cog(merged)
            out_path = Path(f"{args.summary_output_prefix}_cog.{ext}")
            cog_df.to_csv(out_path, sep=args.sep, index=False)
            print(f"COG category summary: {len(cog_df)} samples -> {out_path}")

        if "kegg" in args.summary:
            kegg_df = summarize_kegg(merged)
            out_path = Path(f"{args.summary_output_prefix}_kegg.{ext}")
            kegg_df.to_csv(out_path, sep=args.sep, index=False)
            print(f"KEGG pathway summary: {len(kegg_df)} samples, "
                  f"{kegg_df.shape[1] - 1} pathways -> {out_path}")

        if "completeness" in args.summary:
            completeness_df = summarize_completeness(merged)
            out_path = Path(f"{args.summary_output_prefix}_completeness.{ext}")
            completeness_df.to_csv(out_path, sep=args.sep, index=False)
            print(f"Annotation completeness summary: {len(completeness_df)} samples -> {out_path}")


if __name__ == "__main__":
    main()

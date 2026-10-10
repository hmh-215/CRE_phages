#!/usr/bin/env python3
"""
detect_coverage_breakpoints.py

Scan `samtools depth` per-base output for abrupt coverage anomalies:

  1. Sustained coverage DROPS  -> possible assembly break / mis-join / low-coverage region.
     Reference used: a LOCAL rolling mean depth (window = --window), so drops are judged
     against the coverage immediately around them.

  2. Sustained coverage JUMPS to a high plateau -> possible collapsed repeat / copy-number
     transition boundary. Reference used: the GLOBAL per-contig mean depth taken from
     `samtools coverage` output (coverage_summary.tsv), so jumps are judged against the
     contig's overall expected depth rather than a local window that the jump itself would
     drag upward.

  3. Zero-coverage GAPS (missing positions in the depth file). These only appear if you ran
     `samtools depth` WITHOUT `-a`/`-a -a` (which omit zero-depth positions entirely). Gaps
     at least --drop-size long are reported as drops too, since a silently-missing stretch
     of positions is itself a coverage breakpoint.

Both drop and jump calls require the anomaly to be sustained for a minimum span
(--drop-size / --jump-size) so single-base noise is ignored, as requested.

INPUT ASSUMPTIONS
  --input             : samtools depth -a output, tab-separated, NO header:
                            <contig>\t<pos>\t<depth>
                        Rows for a given contig are expected to be contiguous/sorted
                        (this is how `samtools depth` naturally emits them).
  --coverage-summary  : samtools coverage output (tsv, has a header row starting with
                        '#rname' and a 'meandepth' column). Optional -- if omitted, the
                        per-contig median depth computed from --input is used as the
                        global baseline instead (a robust but slightly weaker fallback).

OUTPUT
  A .tsv or .csv (auto-detected from --output extension) with one row per flagged region:
    contig, start, end, length_bp, type, mean_depth, baseline_depth, baseline_type,
    fraction_or_fold, note

USAGE EXAMPLE
  python detect_coverage_breakpoints.py \\
      --input depth_per_base.tsv \\
      --coverage-summary coverage_summary.tsv \\
      --window 500 \\
      --cov-cutoff 0.2 --drop-size 50 \\
      --jump-cutoff 3.0 --jump-size 50 \\
      --output breakpoints_report.tsv
"""

import argparse
import sys
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Detect abrupt coverage drops/jumps in samtools depth per-base output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True,
                   help="samtools depth per-base tsv (contig<TAB>pos<TAB>depth), no header.")
    p.add_argument("--coverage-summary", default=None,
                   help="samtools coverage tsv, used as the global per-contig baseline "
                        "(column 'meandepth') for jump detection. If omitted, per-contig "
                        "median depth from --input is used instead.")
    p.add_argument("--window", type=int, default=500,
                   help="Rolling window size in bp used to compute the LOCAL baseline mean "
                        "depth for drop detection.")
    p.add_argument("--cov-cutoff", type=float, default=0.2,
                   help="Fraction of the local window baseline depth at/below which a base "
                        "is considered 'low' (e.g. 0.2 = 20%%).")
    p.add_argument("--drop-size", type=int, default=50,
                   help="Minimum contiguous span (bp) of 'low' bases required to report a "
                        "drop. Shorter dips are ignored as noise.")
    p.add_argument("--jump-cutoff", type=float, default=3.0,
                   help="Fold-increase over the GLOBAL contig baseline depth at/above which "
                        "a base is considered part of a high plateau (e.g. 3.0 = 300%%).")
    p.add_argument("--jump-size", type=int, default=50,
                   help="Minimum contiguous span (bp) of 'high' bases required to report a "
                        "jump/plateau. Shorter spikes are ignored as noise.")
    p.add_argument("--min-baseline-depth", type=float, default=5.0,
                   help="Contigs whose global baseline depth is below this are skipped "
                        "entirely (too low-coverage to call breakpoints reliably).")
    p.add_argument("--merge-gap", type=int, default=10,
                   help="Merge two flagged regions of the SAME type on the SAME contig if "
                        "they are separated by less than this many bp (closes small noisy "
                        "gaps inside one real event).")
    p.add_argument("--output", default="breakpoints_report.tsv",
                   help="Output report path. Extension .csv -> comma-separated, otherwise "
                        "tab-separated.")
    return p.parse_args()


def find_runs(mask):
    """Vectorized contiguous-True-run finder. Returns list of (start_idx, end_idx) inclusive."""
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0:
        return []
    d = np.diff(mask.astype(np.int8))
    starts = np.where(d == 1)[0] + 1
    ends = np.where(d == -1)[0]
    if mask[0]:
        starts = np.insert(starts, 0, 0)
    if mask[-1]:
        ends = np.append(ends, mask.size - 1)
    return list(zip(starts.tolist(), ends.tolist()))


def merge_close_runs(runs, pos_arr, merge_gap):
    """Merge runs (index pairs) whose genomic gap between them is < merge_gap bp."""
    if not runs:
        return runs
    runs = sorted(runs)
    merged = [runs[0]]
    for s, e in runs[1:]:
        ps, pe = merged[-1]
        gap = pos_arr[s] - pos_arr[pe] - 1
        if gap < merge_gap:
            merged[-1] = (ps, e)
        else:
            merged.append((s, e))
    return merged


def load_coverage_summary(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#").strip() for c in df.columns]
    missing = {"rname", "meandepth"} - set(df.columns)
    if missing:
        raise ValueError(
            f"--coverage-summary is missing expected column(s) {missing}. "
            "Is this really the output of `samtools coverage`?"
        )
    return dict(zip(df["rname"], df["meandepth"]))


def process_contig(contig, pos_arr, depth_arr, global_baseline, baseline_type, args):
    """Return a list of result dict rows for one contig."""
    rows = []
    n = len(pos_arr)
    if n == 0:
        return rows

    if global_baseline is None:
        global_baseline = float(np.median(depth_arr))
        baseline_type = "median(--input, fallback)"

    if global_baseline < args.min_baseline_depth:
        print(f"  [skip] {contig}: global baseline depth {global_baseline:.2f} < "
              f"--min-baseline-depth ({args.min_baseline_depth}); skipping.", file=sys.stderr)
        return rows

    # ---- 0. Zero-coverage gaps (missing positions => implied zero depth) -----------------
    if n > 1:
        gaps = pos_arr[1:] - pos_arr[:-1] - 1  # missing bp between consecutive rows
        gap_idx = np.where(gaps >= args.drop_size)[0]
        for i in gap_idx:
            gap_start = int(pos_arr[i]) + 1
            gap_end = int(pos_arr[i + 1]) - 1
            length = gap_end - gap_start + 1
            rows.append(dict(
                contig=contig, start=gap_start, end=gap_end, length_bp=length,
                type="zero_coverage_gap",
                mean_depth=0.0, baseline_depth=round(global_baseline, 3),
                baseline_type=baseline_type, fraction_or_fold=0.0,
                note=(f"No depth records for {length} bp (implied zero coverage). "
                      f"If unexpected, re-run `samtools depth -a` to include zero-depth "
                      f"positions explicitly."),
            ))

    # ---- 1. Local rolling mean (reference for drop detection) ----------------------------
    local_mean = pd.Series(depth_arr).rolling(
        window=args.window, min_periods=1, center=True
    ).mean().to_numpy()

    drop_threshold = args.cov_cutoff * local_mean
    drop_mask = depth_arr <= drop_threshold
    drop_runs = find_runs(drop_mask)
    drop_runs = [(s, e) for s, e in drop_runs if (pos_arr[e] - pos_arr[s] + 1) >= args.drop_size]
    drop_runs = merge_close_runs(drop_runs, pos_arr, args.merge_gap)

    for s, e in drop_runs:
        length = int(pos_arr[e] - pos_arr[s] + 1)
        region_mean = float(depth_arr[s:e + 1].mean())
        local_ref = float(local_mean[s:e + 1].mean())
        frac = region_mean / local_ref if local_ref > 0 else float("nan")
        rows.append(dict(
            contig=contig, start=int(pos_arr[s]), end=int(pos_arr[e]), length_bp=length,
            type="coverage_drop",
            mean_depth=round(region_mean, 3), baseline_depth=round(local_ref, 3),
            baseline_type=f"local rolling mean (window={args.window}bp)",
            fraction_or_fold=round(frac, 3),
            note=(f"Depth drops to ~{frac * 100:.1f}% of local window baseline "
                  f"({local_ref:.1f}x) for {length} bp (>= --drop-size {args.drop_size})."),
        ))

    # ---- 2. Global baseline (reference for jump/plateau detection) -----------------------
    jump_threshold = args.jump_cutoff * global_baseline
    jump_mask = depth_arr >= jump_threshold
    jump_runs = find_runs(jump_mask)
    jump_runs = [(s, e) for s, e in jump_runs if (pos_arr[e] - pos_arr[s] + 1) >= args.jump_size]
    jump_runs = merge_close_runs(jump_runs, pos_arr, args.merge_gap)

    for s, e in jump_runs:
        length = int(pos_arr[e] - pos_arr[s] + 1)
        region_mean = float(depth_arr[s:e + 1].mean())
        fold = region_mean / global_baseline if global_baseline > 0 else float("nan")
        rows.append(dict(
            contig=contig, start=int(pos_arr[s]), end=int(pos_arr[e]), length_bp=length,
            type="coverage_jump_plateau",
            mean_depth=round(region_mean, 3), baseline_depth=round(global_baseline, 3),
            baseline_type=baseline_type,
            fraction_or_fold=round(fold, 3),
            note=(f"Depth jumps to ~{fold:.2f}x the contig-wide baseline depth "
                  f"({global_baseline:.1f}x) for {length} bp (>= --jump-size {args.jump_size}). "
                  f"Possible collapsed repeat / copy-number transition boundary."),
        ))

    return rows


def main():
    args = parse_args()

    baseline_map = {}
    if args.coverage_summary:
        print(f"Loading global baselines from {args.coverage_summary} ...", file=sys.stderr)
        baseline_map = load_coverage_summary(args.coverage_summary)
        baseline_type = f"contig meandepth ({args.coverage_summary})"
    else:
        print("No --coverage-summary given; falling back to per-contig median depth "
              "computed from --input for jump detection.", file=sys.stderr)
        baseline_type = "median(--input, fallback)"

    print(f"Reading depth file {args.input} ... (this can take a while for large files)",
          file=sys.stderr)
    df = pd.read_csv(
        args.input, sep="\t", header=None, names=["contig", "pos", "depth"],
        dtype={"pos": "int64", "depth": "int64"},
    )
    print(f"Loaded {len(df):,} rows across {df['contig'].nunique()} contig(s).", file=sys.stderr)

    all_rows = []
    for contig, sub in df.groupby("contig", sort=False):
        pos_arr = sub["pos"].to_numpy()
        depth_arr = sub["depth"].to_numpy()
        order = np.argsort(pos_arr, kind="mergesort")
        pos_arr, depth_arr = pos_arr[order], depth_arr[order]

        global_baseline = baseline_map.get(contig) if baseline_map else None
        if baseline_map and contig not in baseline_map:
            print(f"  [warn] {contig}: not found in --coverage-summary; falling back to "
                  f"per-contig median depth.", file=sys.stderr)

        print(f"Scanning {contig} ({len(pos_arr):,} positions) ...", file=sys.stderr)
        rows = process_contig(
            contig, pos_arr, depth_arr, global_baseline,
            baseline_type if global_baseline is not None else "median(--input, fallback)",
            args,
        )
        all_rows.extend(rows)

    result = pd.DataFrame(all_rows, columns=[
        "contig", "start", "end", "length_bp", "type",
        "mean_depth", "baseline_depth", "baseline_type", "fraction_or_fold", "note",
    ])
    if not result.empty:
        result = result.sort_values(["contig", "start"]).reset_index(drop=True)

    sep = "," if args.output.lower().endswith(".csv") else "\t"
    result.to_csv(args.output, sep=sep, index=False)

    n_drop = (result["type"] == "coverage_drop").sum() if not result.empty else 0
    n_jump = (result["type"] == "coverage_jump_plateau").sum() if not result.empty else 0
    n_gap = (result["type"] == "zero_coverage_gap").sum() if not result.empty else 0
    print(f"\nDone. {len(result)} region(s) flagged "
          f"({n_drop} drop, {n_jump} jump/plateau, {n_gap} zero-coverage gap).",
          file=sys.stderr)
    print(f"Report written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

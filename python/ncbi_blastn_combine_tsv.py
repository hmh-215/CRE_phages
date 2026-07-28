# -*- coding: utf-8 -*-
"""
combine_blast_tsv.py
====================
Combines all individual BLAST TSV output files from a directory
into a single merged TSV file.

Usage
-----
    python combine_blast_tsv.py \
        --input  /path/to/blast_results/ \
        --output /path/to/blast_combined.tsv

The header is written once from the first file found.
All subsequent files have their header line skipped.
Files are merged in alphabetical order by filename.
"""

from __future__ import print_function

import argparse
import glob
import os
import re
import sys


def combine(input_dir, output_file):

    # Natural sort key: splits filename into text and integer chunks so that
    # NODE_2 < NODE_10 < NODE_100 instead of lexicographic NODE_100 < NODE_2
    def natural_key(path):
        parts = re.split(r"(\d+)", os.path.basename(path))
        return [int(p) if p.isdigit() else p.lower() for p in parts]

    # Find all TSV files, excluding any existing summary or combined files
    pattern = os.path.join(input_dir, "*.tsv")
    all_files = sorted(glob.glob(pattern), key=natural_key)

    # Skip summary logs and any previously combined output
    tsv_files = [
        f for f in all_files
        if not os.path.basename(f).startswith("blast_summary")
        and os.path.basename(f) != os.path.basename(output_file)
    ]

    if not tsv_files:
        print("[ERROR] No TSV files found in: {0}".format(input_dir))
        sys.exit(1)

    print("[INFO] Found {0} TSV file(s) to combine.".format(len(tsv_files)))

    rows_written = 0
    header_written = False

    with open(output_file, "w") as out:
        for i, tsv_path in enumerate(tsv_files, 1):
            fname = os.path.basename(tsv_path)
            try:
                with open(tsv_path, "r") as f:
                    lines = f.readlines()

                if not lines:
                    print("  [{0:>3}] SKIP (empty): {1}".format(i, fname))
                    continue

                # Write header only once, from the first non-empty file
                if not header_written:
                    out.write(lines[0])
                    header_written = True

                # Write data rows (skip header line of each file)
                data_lines = [l for l in lines[1:] if l.strip()]
                for line in data_lines:
                    out.write(line)

                rows_written += len(data_lines)
                print("  [{0:>3}] {1:>4} row(s)  {2}".format(
                    i, len(data_lines), fname))

            except Exception as exc:
                print("  [{0:>3}] ERROR reading {1}: {2}".format(i, fname, exc))

    print("\n[DONE] {0} total data rows written to: {1}".format(
        rows_written, output_file))


def main():
    parser = argparse.ArgumentParser(
        description="Combine all BLAST TSV result files into one file."
    )
    parser.add_argument("--input",  required=True,
                        help="Directory containing individual TSV files")
    parser.add_argument("--output", required=True,
                        help="Path for the combined output TSV file")
    args = parser.parse_args()

    if not os.path.isdir(args.input):
        raise IOError("Input directory not found: {0}".format(args.input))

    # Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    combine(args.input, args.output)


if __name__ == "__main__":
    main()

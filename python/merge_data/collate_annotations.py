#!/usr/bin/env python3
"""
collate_annotations.py
----------------------
Collate Panaroo pangenome clusters, Prokka/Bakta GFF coordinates and
eggNOG-mapper functional annotations into one per-gene annotation table,
plus optional cluster and eggNOG (COG / KEGG / completeness) summaries.

BACKGROUND
----------
Locus tags (e.g. GKFICNOF_04064) are unique per isolate, gene symbols are
inconsistent between tools, and Panaroo stores the locus tags of each ortholog
cluster inside per-sample cells. This script joins all three layers on the key
(Sample, Locus_Tag). It replaces merge_eggnog.py + merge_pangenome_annotations.py:
eggNOG files are read directly from a per-sample manifest, so the ~250 MB
concatenated eggNOG table is no longer needed (it can still be written with
--eggnog-full-output), and one sample-name normaliser is applied to eggNOG,
GFF and Panaroo names. That removes the former failure in which eggNOG
samples were labelled '<sample>.prokka' and never matched Panaroo's '<sample>'.

KEY CONCEPTS
------------
1. Modes
     Panaroo mode  (--panaroo given): one row per locus tag listed in the
                   Panaroo matrix, with cluster ID, compartment and frequency,
                   GFF coordinates (Prokka GFFs used for Panaroo) and eggNOG.
     GFF mode      (no --panaroo): one row per CDS locus tag in the GFFs
                   (e.g. Bakta GFF3 + Bakta eggNOG).
2. Compartments (cohort frequency f): Core f >= 0.99, SoftCore >= 0.95,
   Shell >= 0.15, Cloud < 0.15 (cutoffs adjustable). With --samples the
   frequency is computed over cohort genomes only.
3. eggNOG join: exact (Sample, Locus_Tag) first; if the sample name is not
   found, a unique locus-tag-only match is accepted (Prokka/Bakta tags carry a
   per-genome prefix). The summary reports the match rate; a low rate almost
   always means eggNOG was run on a different annotation (Prokka vs Bakta)
   than the GFFs / Panaroo input.
4. Panaroo cells may hold several tags (separated by ';', ',' or TAB) and
   'ref~~~tag' prefixes; every tag is expanded. Trailing _stop / _len / _pseudo
   flags are tolerated when matching GFF and eggNOG.

Directory conventions (matching Ecoli_merge_script_TN.bash)
-----------------------------------------------------------
  cohort list          : {PATHS_DIR}/samples_TN.txt
  Panaroo matrix       : {PANAROO_RES_DIR}/gene_presence_absence.csv
  Prokka GFFs          : {PANAROO_GFF_DIR}/<sample>.gff
  Bakta GFF manifest   : {PATHS_DIR}/paths_bakta_gff_TN.txt
  eggNOG manifests     : {PATHS_DIR}/paths_eggnog_{prokka,bakta}_TN.txt
  per-gene tables      : {METADATA_OUT_DIR}/*_annotated_genes_TN.tsv (+ .summary.txt)
  eggNOG summaries     : {INDIVIDUAL_OUT_DIR}/eggnog_{prokka,bakta}_{cog,kegg,completeness}.tsv

Usage examples
--------------
# Example 1: Panaroo + Prokka GFF + Prokka eggNOG, with cluster and eggNOG summaries
python collate_annotations.py \\
    --output          merge_metadata/panaroo_prokka_eggnog_annotated_genes_TN.tsv \\
    --panaroo         panaroo_Ecoli_TN/results/gene_presence_absence.csv \\
    --gff-dir         panaroo_Ecoli_TN/gff \\
    --eggnog-manifest paths/paths_eggnog_prokka_TN.txt \\
    --samples         paths/samples_TN.txt \\
    --cluster-summary merge_metadata/panaroo_clusters_summary_TN.tsv \\
    --eggnog-summary  cog kegg completeness \\
    --eggnog-summary-prefix merge_individual_outputs/eggnog_prokka

# Example 2: Bakta GFF3 + Bakta eggNOG (no pangenome)
python collate_annotations.py \\
    --output          merge_metadata/bakta_eggnog_annotated_genes_TN.tsv \\
    --gff-manifest    paths/paths_bakta_gff_TN.txt \\
    --eggnog-manifest paths/paths_eggnog_bakta_TN.txt \\
    --samples         paths/samples_TN.txt
"""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import unquote


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

PANAROO_META_COLS = {
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC",
    "Min group size nuc", "Max group size nuc", "Avg group size nuc",
}

# eggNOG fields kept per gene: (output column, eggNOG column).
EGGNOG_FIELDS = [
    ("eggNOG_Preferred_Name", "Preferred_name"),
    ("eggNOG_Description", "Description"),
    ("COG_Category", "COG_category"),
    ("GO_Terms", "GOs"),
    ("KEGG_ko", "KEGG_ko"),
    ("KEGG_Pathway", "KEGG_Pathway"),
]

# Columns checked for "annotated" in the completeness summary.
COMPLETENESS_FIELDS = [
    ("seed_ortholog", "Any_hit"),
    ("COG_category", "COG_category"),
    ("GOs", "GO_terms"),
    ("KEGG_ko", "KEGG_ortholog"),
    ("PFAMs", "PFAM"),
]

TAG_FLAG_RE = re.compile(r"_(stop|len|pseudo)$")


def normalise(text: str) -> str:
    """Lowercase and strip hyphens, underscores, dots and spaces for fuzzy matching."""
    return re.sub(r"[-_ .]", "", str(text).lower())


def clean_sample_name(raw: str) -> str:
    """
    Reduce a file name, path or tool label to a bare sample ID.

    e.g. '132895-LJF30331.prokka.eggnog.emapper.annotations' -> '132895-LJF30331'.

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


def clean_locus_tag(raw: str) -> str:
    """Strip Panaroo 'ref~~~' prefixes and whitespace from a locus tag."""
    tag = raw.strip()
    if "~~~" in tag:
        tag = tag.split("~~~")[-1].strip()
    return tag


def split_cell_tags(cell: str) -> list[str]:
    """Expand a Panaroo sample cell into individual cleaned locus tags."""
    if not cell or cell.strip().lower() in ("", "-", "na", "n/a"):
        return []
    parts = re.split(r"[;,\t]", cell)
    return [t for t in (clean_locus_tag(p) for p in parts) if t]


def tag_variants(tag: str) -> list[str]:
    """Return the tag plus its form without a trailing _stop/_len/_pseudo flag."""
    stripped = TAG_FLAG_RE.sub("", tag)
    return [tag] if stripped == tag else [tag, stripped]


def classify_compartment(freq: float, core: float, soft: float, shell: float) -> str:
    """Assign Core / SoftCore / Shell / Cloud from cohort frequency (0-1)."""
    if freq >= core:
        return "Core"
    if freq >= soft:
        return "SoftCore"
    if freq >= shell:
        return "Shell"
    return "Cloud"


def resolve_annotations_file(path: Path) -> Path | None:
    """Return path if it is a file, else the first *.emapper.annotations inside it."""
    if path.is_file():
        return path
    if path.is_dir():
        hits = sorted(path.glob("*.emapper.annotations")) or sorted(path.rglob("*.emapper.annotations"))
        return hits[0] if hits else None
    return None


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def read_manifest(manifest_path: Path) -> list[tuple[str, Path]]:
    """
    Parse a manifest of 'sample<TAB>path' or bare-path lines (TAB-only splitting).

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
    """Read the cohort list (first TAB field per non-comment line)."""
    if path is None:
        return []
    if not path.exists():
        sys.exit(f"ERROR: --samples file not found: {path}")
    out: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line and not line.startswith("#"):
                s = clean_sample_name(line.split("\t")[0])
                if s not in out:
                    out.append(s)
    return out


def parse_gff(path: Path) -> dict[str, dict[str, str]]:
    """
    Map locus_tag -> coordinates and annotation from a Prokka GFF or Bakta GFF3.

    Only CDS/gene features are used; parsing stops at '##FASTA'. CDS
    coordinates override gene coordinates for the same tag.

    Returns:
        {locus_tag: {"contig", "start", "end", "strand", "gene", "product"}}.
    """
    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith("##FASTA"):
                break
            if raw.startswith("#") or not raw.strip():
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("CDS", "gene"):
                continue
            attrs: dict[str, str] = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k.strip()] = unquote(v.strip())
            tag = clean_locus_tag(attrs.get("locus_tag") or attrs.get("ID") or "")
            if not tag:
                continue
            if tag not in out or parts[2] == "CDS":
                out[tag] = {
                    "contig": parts[0], "start": parts[3], "end": parts[4], "strand": parts[6],
                    "gene": attrs.get("gene") or attrs.get("Name") or "",
                    "product": attrs.get("product") or "",
                }
    return out


def iter_eggnog_file(path: Path):
    """
    Yield (header, row_fields) from an eggNOG annotations file or merged table.

    '##' comment lines are skipped; the header may start with '#query' or 'query'
    (merged tables produced by merge_eggnog.py start with 'Sample').
    """
    header: list[str] | None = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith("##") or not raw.strip():
                continue
            line = raw.rstrip("\n")
            if header is None:
                header = line.lstrip("#").split("\t")
                continue
            fields = line.split("\t")
            if len(fields) == len(header):
                yield header, fields


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

class EggnogIndex:
    """
    Locus-tag index of eggNOG annotations with per-sample summary accumulators.

    Stores one tuple of the EGGNOG_FIELDS values per (sample, tag); keyed by tag
    first so a unique tag-only match can rescue residual sample-name mismatches.
    """

    def __init__(self, summaries: set[str]):
        """Initialise empty index and the requested summary accumulators."""
        self.by_tag: dict[str, list[tuple[str, tuple[str, ...]]]] = {}
        self.samples: list[str] = []
        self.summaries = summaries
        self.cog: dict[str, Counter] = {}
        self.kegg: dict[str, Counter] = {}
        self.complete: dict[str, Counter] = {}
        self.total: Counter = Counter()

    def add(self, sample: str, rec: dict[str, str]) -> None:
        """Index one eggNOG row and update summary counters."""
        tag = clean_locus_tag(rec.get("query", "") or rec.get("query_name", "") or "")
        if not tag:
            return
        if sample not in self.total:
            self.samples.append(sample)
        self.total[sample] += 1
        values = tuple((rec.get(src, "") or "").replace("-", "") if (rec.get(src, "") or "") == "-"
                       else (rec.get(src, "") or "") for _, src in EGGNOG_FIELDS)
        self.by_tag.setdefault(tag, []).append((sample, values))
        if "cog" in self.summaries:
            cog = (rec.get("COG_category", "") or "").strip()
            c = self.cog.setdefault(sample, Counter())
            if cog in ("", "-"):
                c["Unclassified"] += 1
            else:
                for letter in cog:
                    c[letter] += 1
        if "kegg" in self.summaries:
            pw = (rec.get("KEGG_Pathway", "") or "").strip()
            k = self.kegg.setdefault(sample, Counter())
            if pw in ("", "-"):
                k["Unassigned"] += 1
            else:
                for p in pw.split(","):
                    if p.strip():
                        k[p.strip()] += 1
        if "completeness" in self.summaries:
            cc = self.complete.setdefault(sample, Counter())
            for col, label in COMPLETENESS_FIELDS:
                if col in rec and (rec.get(col, "") or "-") != "-":
                    cc[label] += 1

    def lookup(self, sample: str, tag: str) -> tuple[tuple[str, ...] | None, str]:
        """
        Find eggNOG values for (sample, tag).

        Returns:
            (values or None, match kind: "exact" | "tag_only" | "none").
        """
        s_norm = normalise(sample)
        for t in tag_variants(tag):
            hits = self.by_tag.get(t)
            if not hits:
                continue
            for hs, vals in hits:
                if hs == sample or normalise(hs) == s_norm:
                    return vals, "exact"
            if len(hits) == 1:
                return hits[0][1], "tag_only"
        return None, "none"


def load_eggnog(
    manifest: Path | None,
    table: Path | None,
    cohort: set[str],
    summaries: set[str],
    full_output: Path | None,
) -> EggnogIndex:
    """
    Read eggNOG annotations from a per-sample manifest and/or a merged table.

    Args:
        manifest: 'sample<TAB>path' manifest of .emapper.annotations files/dirs.
        table: Legacy merged table with a 'Sample' column (merge_eggnog.py output).
        cohort: Samples to keep (empty = all).
        summaries: Summary types to accumulate.
        full_output: If set, stream every kept row to this concatenated TSV.

    Returns:
        Populated EggnogIndex.
    """
    idx = EggnogIndex(summaries)
    writer = None
    out_fh = None
    sources: list[tuple[str | None, Path]] = []
    if manifest:
        for sample, path in read_manifest(manifest):
            f = resolve_annotations_file(path)
            if f is None:
                print(f"  [{sample}] WARNING: no .emapper.annotations at {path} - skipping", file=sys.stderr)
                continue
            sources.append((sample, f))
    if table:
        if table.is_file():
            sources.append((None, table))
        else:
            print(f"WARNING: eggNOG table not found: {table}", file=sys.stderr)

    try:
        for sample, f in sources:
            for header, fields in iter_eggnog_file(f):
                rec = dict(zip(header, fields))
                s = sample or clean_sample_name(rec.get("Sample", "") or "")
                if cohort and s not in cohort:
                    continue
                idx.add(s, rec)
                if full_output:
                    if writer is None:
                        full_output.parent.mkdir(parents=True, exist_ok=True)
                        out_fh = open(full_output, "w", newline="", encoding="utf-8")
                        cols = ["Sample"] + [h for h in header if h != "Sample"]
                        writer = csv.DictWriter(out_fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
                        writer.writeheader()
                    rec["Sample"] = s
                    writer.writerow(rec)
    finally:
        if out_fh:
            out_fh.close()
            print(f"  Saved: {full_output} (concatenated eggNOG rows)")
    print(f"Loaded eggNOG annotations: {sum(idx.total.values())} proteins from {len(idx.samples)} samples.")
    return idx


def collect_gffs(gff_dir: Path | None, gff_manifest: Path | None, cohort: set[str]) -> dict[str, Path]:
    """
    Collect one GFF path per sample from a directory and/or manifest.

    Returns:
        {sample: gff_path}; manifest entries take precedence over directory hits.
    """
    out: dict[str, Path] = {}
    if gff_dir and gff_dir.is_dir():
        for pat in ("*.gff", "*.gff3", "*/*.gff", "*/*.gff3"):
            for p in sorted(gff_dir.glob(pat)):
                out.setdefault(sample_from_path(p), p)
    if gff_manifest:
        for s, p in read_manifest(gff_manifest):
            if p.is_file():
                out[s] = p
            else:
                print(f"  [{s}] WARNING: GFF not found: {p}", file=sys.stderr)
    if cohort:
        out = {s: p for s, p in out.items() if s in cohort}
    return out


def egg_columns(vals: tuple[str, ...] | None) -> dict[str, str]:
    """Expand an eggNOG value tuple into output columns (blank when unmatched)."""
    if vals is None:
        return {col: "" for col, _ in EGGNOG_FIELDS}
    return {col: v for (col, _), v in zip(EGGNOG_FIELDS, vals)}


def collate_panaroo(
    panaroo: Path,
    gffs: dict[str, Path],
    egg: EggnogIndex | None,
    cohort: list[str],
    cutoffs: tuple[float, float, float],
    stats: Counter,
) -> tuple[list[dict], list[dict], int]:
    """
    Join Panaroo clusters with GFF coordinates and eggNOG annotations.

    Returns:
        (per-gene rows, cluster summary rows, number of genomes used).
    """
    with open(panaroo, encoding="utf-8", errors="replace") as fh:
        delim = "\t" if "\t" in fh.readline() else ","
    gff_cache = {s: parse_gff(p) for s, p in gffs.items()}
    gff_norm = {normalise(s): v for s, v in gff_cache.items()}
    genes: list[dict] = []
    clusters: list[dict] = []
    with open(panaroo, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        col_to_sample = {c: clean_sample_name(c) for c in (reader.fieldnames or []) if c not in PANAROO_META_COLS}
        if cohort:
            keep = set(cohort)
            col_to_sample = {c: s for c, s in col_to_sample.items() if s in keep}
        n = len(col_to_sample)
        if n == 0:
            sys.exit("ERROR: no usable sample columns in the Panaroo matrix (check --samples).")
        for row in reader:
            per_sample = {s: split_cell_tags(row.get(c, "") or "") for c, s in col_to_sample.items()}
            per_sample = {s: t for s, t in per_sample.items() if t}
            if not per_sample:
                continue
            freq = len(per_sample) / n
            comp = classify_compartment(freq, *cutoffs)
            cluster = row.get("Gene", "") or ""
            clusters.append({
                "Panaroo_Cluster": cluster, "Pangenome_Class": comp,
                "Frequency_Pct": round(freq * 100, 2), "Isolates_Count": len(per_sample),
                "Total_Genomes": n, "Consensus_Gene": row.get("Non-unique Gene name", "") or "",
                "Consensus_Annotation": row.get("Annotation", "") or "",
            })
            for s, tags in per_sample.items():
                coords = gff_cache.get(s) or gff_norm.get(normalise(s)) or {}
                for tag in tags:
                    g = next((coords[t] for t in tag_variants(tag) if t in coords), {})
                    stats["gff_hit" if g else "gff_miss"] += 1
                    vals, kind = egg.lookup(s, tag) if egg else (None, "none")
                    stats[f"egg_{kind}"] += 1
                    genes.append({
                        "Sample": s, "Locus_Tag": tag, "Panaroo_Cluster": cluster,
                        "Pangenome_Class": comp, "Cluster_Freq_Pct": round(freq * 100, 2),
                        "Contig": g.get("contig", ""), "Start": g.get("start", ""),
                        "End": g.get("end", ""), "Strand": g.get("strand", ""),
                        "Annotation_Gene": g.get("gene", ""), "Annotation_Product": g.get("product", ""),
                        **egg_columns(vals),
                    })
    return genes, clusters, n


def collate_gff_only(
    gffs: dict[str, Path],
    egg: EggnogIndex | None,
    stats: Counter,
) -> list[dict]:
    """
    Join GFF CDS coordinates with eggNOG annotations (no pangenome).

    Returns:
        Per-gene rows.
    """
    genes: list[dict] = []
    for s in sorted(gffs):
        for tag, g in parse_gff(gffs[s]).items():
            vals, kind = egg.lookup(s, tag) if egg else (None, "none")
            stats[f"egg_{kind}"] += 1
            genes.append({
                "Sample": s, "Locus_Tag": tag,
                "Contig": g["contig"], "Start": g["start"], "End": g["end"], "Strand": g["strand"],
                "Annotation_Gene": g["gene"], "Annotation_Product": g["product"],
                **egg_columns(vals),
            })
    return genes


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

CLUSTER_COLUMNS = ["Panaroo_Cluster", "Pangenome_Class", "Frequency_Pct", "Isolates_Count",
                   "Total_Genomes", "Consensus_Gene", "Consensus_Annotation"]


def gene_columns(has_panaroo: bool, has_gff: bool, has_eggnog: bool) -> list[str]:
    """Build per-gene output columns, pruning blocks whose inputs were omitted."""
    cols = ["Sample", "Locus_Tag"]
    if has_panaroo:
        cols += ["Panaroo_Cluster", "Pangenome_Class", "Cluster_Freq_Pct"]
    if has_gff:
        cols += ["Contig", "Start", "End", "Strand", "Annotation_Gene", "Annotation_Product"]
    if has_eggnog:
        cols += [c for c, _ in EGGNOG_FIELDS]
    return cols


def write_table(rows: list[dict], columns: list[str], outpath: Path, comments: list[str] | None = None) -> None:
    """Write rows to a TSV with optional leading '# ' comment lines."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        for c in comments or []:
            fh.write(f"# {c}\n")
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {outpath} ({len(rows)} rows x {len(columns)} columns)")


def write_counter_matrix(data: dict[str, Counter], last: str, outpath: Path, comment: str) -> None:
    """Write a sample x category count matrix, keeping the catch-all column last."""
    cats = sorted({k for c in data.values() for k in c if k != last})
    if any(last in c for c in data.values()):
        cats.append(last)
    rows = [{"Sample": s, **{k: data[s].get(k, 0) for k in cats}} for s in sorted(data)]
    write_table(rows, ["Sample"] + cats, outpath, [comment])


def write_eggnog_summaries(egg: EggnogIndex, prefix: Path) -> None:
    """Write the requested COG / KEGG / completeness summaries."""
    if "cog" in egg.summaries:
        write_counter_matrix(egg.cog, "Unclassified", Path(f"{prefix}_cog.tsv"),
                             "TOOL: collate_annotations.py (COG category letter counts per sample; multi-letter categories count once per letter)")
    if "kegg" in egg.summaries:
        write_counter_matrix(egg.kegg, "Unassigned", Path(f"{prefix}_kegg.tsv"),
                             "TOOL: collate_annotations.py (genes per KEGG pathway per sample)")
    if "completeness" in egg.summaries:
        rows = []
        for s in sorted(egg.total):
            tot = egg.total[s]
            rec = {"Sample": s, "Total_genes": tot}
            for _, label in COMPLETENESS_FIELDS:
                n = egg.complete.get(s, Counter()).get(label, 0)
                rec[f"{label}_n"] = n
                rec[f"{label}_pct"] = round(100 * n / tot, 2) if tot else 0.0
            rows.append(rec)
        cols = ["Sample", "Total_genes"] + [f"{l}_{k}" for _, l in COMPLETENESS_FIELDS for k in ("n", "pct")]
        write_table(rows, cols, Path(f"{prefix}_completeness.tsv"),
                    ["TOOL: collate_annotations.py (eggNOG annotation completeness per sample)"])


def write_summary(
    genes: list[dict],
    clusters: list[dict],
    stats: Counter,
    has_eggnog: bool,
    missing_gff: list[str],
    outpath: Path,
) -> None:
    """Write and echo the annotation reconciliation audit, including eggNOG match rate."""
    total = len(genes)
    lines = [
        "=============================================================",
        " Pangenome / Annotation Collation Summary",
        "=============================================================",
        f"  Samples                    : {len({g['Sample'] for g in genes})}",
        f"  Locus tags (rows)          : {total}",
    ]
    if clusters:
        cc = Counter(c["Pangenome_Class"] for c in clusters)
        lines.append(f"  Panaroo clusters           : {len(clusters)}")
        for k in ("Core", "SoftCore", "Shell", "Cloud"):
            lines.append(f"    {k:<9}                : {cc.get(k, 0)} ({cc.get(k, 0) * 100 / len(clusters):.1f}%)")
        hit = stats.get("gff_hit", 0)
        lines.append(f"  GFF coordinate match       : {hit}/{total} ({hit * 100 / total:.1f}%)" if total else "")
        if missing_gff:
            lines.append(f"  Samples without GFF        : {len(missing_gff)} {missing_gff}")
    if has_eggnog and total:
        ex, to = stats.get("egg_exact", 0), stats.get("egg_tag_only", 0)
        rate = (ex + to) * 100 / total
        lines += [
            f"  eggNOG matched (exact)     : {ex} ({ex * 100 / total:.1f}%)",
            f"  eggNOG matched (tag only)  : {to} ({to * 100 / total:.1f}%)",
            f"  eggNOG unmatched           : {stats.get('egg_none', 0)} ({100 - rate:.1f}%)",
        ]
        if rate < 50:
            lines.append("  WARNING: <50% eggNOG match - eggNOG was probably run on a different "
                         "annotation (Prokka vs Bakta) than the GFF/Panaroo input.")
    lines.append("=============================================================")
    lines = [ln for ln in lines if ln]
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
            "Collate Panaroo clusters, Prokka/Bakta GFF coordinates and eggNOG-mapper annotations "
            "into a per-gene annotation table with optional cluster and eggNOG summaries."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Required ---
    parser.add_argument("--output", required=True, type=Path, metavar="TSV",
                        help="Destination of the per-gene table; <output>.summary.txt is written alongside "
                             "(orchestrator: ${METADATA_OUT_DIR}/*_annotated_genes_TN.tsv).")

    # --- Inputs ---
    parser.add_argument("--panaroo", type=Path, default=None, metavar="CSV",
                        help="Panaroo gene_presence_absence.csv. If given, runs Panaroo mode (one row per "
                             "clustered locus tag); otherwise GFF mode. Orchestrator: ${panaroo_matrix}.")
    parser.add_argument("--gff-dir", type=Path, default=None, metavar="DIR",
                        help="Directory of per-sample .gff/.gff3 files (e.g. the staged Prokka GFFs used for "
                             "Panaroo). Orchestrator: ${PANAROO_GFF_DIR}.")
    parser.add_argument("--gff-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample GFF paths ('sample<TAB>path'). Overrides --gff-dir hits "
                             "for the same sample. Orchestrator: ${PATHS_DIR}/paths_bakta_gff_TN.txt.")
    parser.add_argument("--eggnog-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample eggNOG .emapper.annotations files or output directories. "
                             "Must come from the SAME annotation as the GFFs (Prokka with Prokka, Bakta with Bakta).")
    parser.add_argument("--eggnog-table", type=Path, default=None, metavar="TSV",
                        help="Legacy merged eggNOG table with a 'Sample' column (merge_eggnog.py --full-output). "
                             "Sample names are normalised, so '<sample>.prokka' labels match.")
    parser.add_argument("--samples", type=Path, default=None, metavar="TXT",
                        help="Cohort list; restricts Panaroo columns, GFFs and eggNOG rows to these isolates and "
                             "computes compartment frequencies over them. Orchestrator: ${PATHS_DIR}/samples_TN.txt.")

    # --- Optional outputs ---
    parser.add_argument("--cluster-summary", type=Path, default=None, metavar="TSV",
                        help="Panaroo mode: write one row per cluster with compartment, frequency and consensus "
                             "name/annotation.")
    parser.add_argument("--eggnog-summary", nargs="+", choices=["cog", "kegg", "completeness"], default=None,
                        help="eggNOG per-sample summaries to write (requires --eggnog-summary-prefix).")
    parser.add_argument("--eggnog-summary-prefix", type=Path, default=None, metavar="PREFIX",
                        help="Prefix for summaries: <PREFIX>_cog.tsv, _kegg.tsv, _completeness.tsv. "
                             "Orchestrator: ${INDIVIDUAL_OUT_DIR}/eggnog_prokka or eggnog_bakta.")
    parser.add_argument("--eggnog-full-output", type=Path, default=None, metavar="TSV",
                        help="Optional: also stream every eggNOG row (all columns + Sample) to this concatenated "
                             "table. Large (~250 MB for 71 genomes); off by default.")

    # --- Cutoffs ---
    parser.add_argument("--core-cutoff", type=float, default=0.99, metavar="FLOAT",
                        help="Cohort frequency for Core genes. Default: 0.99.")
    parser.add_argument("--soft-core-cutoff", type=float, default=0.95, metavar="FLOAT",
                        help="Cohort frequency for Soft-core genes. Default: 0.95.")
    parser.add_argument("--shell-cutoff", type=float, default=0.15, metavar="FLOAT",
                        help="Cohort frequency for Shell genes (below = Cloud). Default: 0.15.")

    return parser.parse_args()


def main():
    """Orchestrate loading, joining and writing."""
    args = parse_args()

    has_eggnog_input = bool(args.eggnog_manifest or args.eggnog_table)
    if args.eggnog_summary and not args.eggnog_summary_prefix:
        sys.exit("ERROR: --eggnog-summary requires --eggnog-summary-prefix.")
    if not args.panaroo and not (args.gff_dir or args.gff_manifest):
        sys.exit("ERROR: provide --panaroo and/or --gff-dir / --gff-manifest.")
    if args.panaroo and not args.panaroo.is_file():
        sys.exit(f"ERROR: Panaroo matrix not found: {args.panaroo}")

    cohort = read_sample_list(args.samples)
    cohort_set = set(cohort)
    mode = "Panaroo" if args.panaroo else "GFF"

    print("\n========================================================")
    print(" collate_annotations.py - per-gene annotation collation")
    print("========================================================")
    print(f" Mode            : {mode}")
    print(f" Panaroo matrix  : {args.panaroo or 'None'}")
    print(f" GFF dir         : {args.gff_dir or 'None'}")
    print(f" GFF manifest    : {args.gff_manifest or 'None'}")
    print(f" eggNOG manifest : {args.eggnog_manifest or 'None'}")
    print(f" eggNOG table    : {args.eggnog_table or 'None'}")
    print(f" Cohort list     : {args.samples or 'None'}{f' ({len(cohort)} samples)' if cohort else ''}")
    print(f" Output          : {args.output}")
    print("========================================================\n")

    egg = None
    if has_eggnog_input:
        egg = load_eggnog(args.eggnog_manifest, args.eggnog_table, cohort_set,
                          set(args.eggnog_summary or []), args.eggnog_full_output)
    gffs = collect_gffs(args.gff_dir, args.gff_manifest, cohort_set)
    print(f"GFF files collected: {len(gffs)}")

    stats: Counter = Counter()
    clusters: list[dict] = []
    if args.panaroo:
        genes, clusters, n_gen = collate_panaroo(
            args.panaroo, gffs, egg, cohort,
            (args.core_cutoff, args.soft_core_cutoff, args.shell_cutoff), stats)
        missing_gff = sorted({g["Sample"] for g in genes} - set(gffs))
    else:
        genes = collate_gff_only(gffs, egg, stats)
        missing_gff = []
        n_gen = len(gffs)
    if not genes:
        sys.exit("ERROR: no gene records produced - check inputs and --samples.")

    has_gff = bool(gffs)
    comments = [
        f"TOOL: collate_annotations.py ({'Panaroo + GFF + eggNOG' if args.panaroo else 'GFF + eggNOG'} per-gene annotations)",
        f"SOURCES: panaroo={args.panaroo.name if args.panaroo else 'omitted'}; "
        f"gff={'included (' + str(len(gffs)) + ' files)' if has_gff else 'omitted'}; "
        f"eggnog={'included' if egg else 'omitted'}",
        "JOIN KEY: (Sample, Locus_Tag); eggNOG fields are blank when no eggNOG hit exists for that locus tag.",
    ]
    if args.panaroo:
        comments.append(f"COMPARTMENTS over {n_gen} genomes: Core >= {args.core_cutoff * 100:.0f}%, "
                        f"SoftCore >= {args.soft_core_cutoff * 100:.0f}%, Shell >= {args.shell_cutoff * 100:.0f}%, "
                        f"Cloud < {args.shell_cutoff * 100:.0f}%")
    write_table(genes, gene_columns(bool(args.panaroo), has_gff, bool(egg)), args.output, comments)

    if args.cluster_summary and clusters:
        write_table(clusters, CLUSTER_COLUMNS, args.cluster_summary, [
            "TOOL: collate_annotations.py (Panaroo cluster summary)",
            comments[-1],
        ])
    if egg and args.eggnog_summary:
        write_eggnog_summaries(egg, args.eggnog_summary_prefix)

    write_summary(genes, clusters, stats, bool(egg), missing_gff, args.output.with_suffix(".summary.txt"))
    print("\nDone.")


if __name__ == "__main__":
    main()

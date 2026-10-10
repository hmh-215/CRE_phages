#!/usr/bin/env python3
"""
collate_genome_overview.py
--------------------------
Collate per-isolate assembly QC (QUAST, CheckM), pangenome composition
(Panaroo), sequence typing (MLST), serotyping (ECTyper or Kleborate),
mobile-genetic-element counts (ISEScan, IntegronFinder, PhiSpy),
plasmid identification (Platon, PlasmidFinder, MOB-suite),
coverage anomaly breakpoints (detect_coverage_breakpoints.py), and
annotation ncRNAs (Prokka/Bakta rRNA, tRNA, repeats) into one master
genome-overview table with optional, user-defined QC flags.

BACKGROUND
----------
Before phage-susceptible and phage-resistant cohorts are compared, every
isolate must be characterised along orthogonal axes: assembly continuity,
marker-gene completeness/contamination, pangenome composition, clonal
background (ST), surface antigens (O:H or K/O loci), mobile-element burden,
and plasmid content.

This script reads raw tool outputs directly, applies one sample-name
normaliser to every source, and restricts the table to an explicit cohort
so batch-wide tables cannot inject extra rows.

KEY CONCEPTS
------------
1. Cohort restriction (--samples): only listed isolates are reported, in list
   order. Isolates found in any source but not in the list are dropped and
   reported in the summary; listed isolates with no data still get a row.
2. QC flags are OFF unless a threshold is given. Each flag adds one column:
     CheckM_Status    PASS  completeness >= --min-completeness and
                            contamination <= --max-contamination
                      WARNING completeness >= 90 and contamination <= 10
                      FAIL  otherwise
     Assembly_Status  PASS <= --max-contigs; WARNING <= max+200; FAIL otherwise
     Size_Status      (--size-status + --species) PASS inside the expected
                      species range; WARNING within +/-10 %; FAIL otherwise
     Pangenome_Status PASS Core_% >= --min-core-pct; WARNING >= value-10;
                      OUTLIER otherwise
     Overall_QC       FAIL if any FAIL/OUTLIER, WARNING if any WARNING, else PASS
   Missing data yields UNKNOWN for that flag.
3. Mobile genetic elements & Plasmids:
     ISEScan        : <genome>.fasta.csv, one row per IS.
     IntegronFinder : <genome>.integrons, one row per element.
     PhiSpy         : prophage_coordinates.tsv / prophage.tsv.
     Platon         : <sample>.tsv (ID, Length, Coverage, Circular, Type).
     PlasmidFinder  : results_tab.tsv or ABRicate plasmidfinder table.
     MOB-suite      : contig_report.txt (mob_recon / mob_typer).
     Breakpoints    : breakpoints_report.tsv (detect_coverage_breakpoints.py).
     Annotations    : Prokka/Bakta summary .txt or GFF (rRNA, tRNA, repeat_region).
4. CheckM_Strain_Het is copied from the CheckM 'Strain heterogeneity' column.
"""

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# File-name tokens stripped (right to left) when inferring a sample name.
SAMPLE_SUFFIX_TOKENS = {
    "tsv", "csv", "txt", "gff", "gff3", "gbk", "fasta", "fna", "fa", "faa",
    "annotations", "emapper", "eggnog", "prokka", "bakta",
    "abricate", "card", "vfdb", "ncbi", "resfinder", "plasmidfinder", "amrfinder",
    "quast", "raw", "filtered", "contigs", "scaffolds",
    "isescan", "integron", "integrons", "summary", "phispy",
    "ectyper", "mlst", "kleborate", "checkm", "platon", "mobsuite", "mobtyper",
    "breakpoints", "chromosome", "plasmid",
}

# Generic per-tool file names that carry no sample ID (use the parent directory).
GENERIC_BASENAMES = {
    "report", "output", "prophage", "prophage_coordinates", "summary", "results",
    "contig_report", "mobtyper_results", "breakpoints_report", "results_tab"
}

SPECIES_SIZE_RANGES = {
    "ecoli": (4_500_000, 5_700_000),
    "klebsiella": (5_000_000, 6_000_000),
    "pseudomonas": (6_000_000, 7_300_000),
}

PANAROO_META_COLS = {
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
    "Avg sequences per isolate", "Genome Fragment", "Order within Fragment",
    "Accessory Fragment", "Accessory Order with Fragment", "QC",
    "Min group size nuc", "Max group size nuc", "Avg group size nuc",
}

MISSING_TOKENS = {"", "-", "na", "n/a", "none", "null", "?", "nd"}


def normalise(text: str) -> str:
    """Lowercase and strip hyphens, underscores, dots and spaces for fuzzy matching."""
    return re.sub(r"[-_ .]", "", str(text).lower())


def clean_sample_name(raw: str) -> str:
    """Reduce a file name, path or tool label to a bare sample ID."""
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
    """Infer a sample ID from a path, falling back to parent directory for generic names."""
    name = clean_sample_name(path.name)
    if name.lower() in GENERIC_BASENAMES:
        name = clean_sample_name(path.parent.name)
    return name


def to_float(val: str | None) -> float | None:
    """Convert a string to float, returning None for missing/unparseable values."""
    if val is None:
        return None
    s = str(val).strip().replace("%", "")
    if s.lower() in MISSING_TOKENS:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def to_int(val: str | None) -> int | None:
    """Convert a string to int (via float), returning None for missing values."""
    f = to_float(val)
    return int(round(f)) if f is not None else None


def pick(row: dict, exact: tuple[str, ...], suffixes: tuple[str, ...] = ()) -> str:
    """Return first non-empty value from exact column names, then from suffix matches."""
    for k in exact:
        v = (row.get(k, "") or "").strip()
        if v:
            return v
    for col, v in row.items():
        if col and v and any(col.endswith(sfx) for sfx in suffixes):
            if str(v).strip():
                return str(v).strip()
    return ""


def format_st(st: str) -> str:
    """Prefix numeric STs with 'ST'; keep '-', novel or partial calls verbatim."""
    st = (st or "").strip()
    if not st or st.lower() in MISSING_TOKENS:
        return "-"
    return f"ST{st}" if st.isdigit() else st


def classify_compartment(freq: float, core: float, soft: float, shell: float) -> str:
    """Assign Core / SoftCore / Shell / Cloud from cohort frequency (0-1)."""
    if freq >= core:
        return "Core"
    if freq >= soft:
        return "SoftCore"
    if freq >= shell:
        return "Shell"
    return "Cloud"


def first_file(base: Path, patterns: tuple[str, ...]) -> Path | None:
    """Return base if it is a file, else the first recursive glob hit under base."""
    if base.is_file():
        return base
    if base.is_dir():
        for pat in patterns:
            hits = sorted(base.rglob(pat))
            if hits:
                return hits[0]
    return None


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def read_manifest(manifest_path: Path) -> list[tuple[str, Path]]:
    """Parse a manifest of 'sample<TAB>path' or bare-path lines."""
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


def expand_table_inputs(paths: list[Path]) -> list[tuple[str | None, Path]]:
    """Expand typing inputs that may be tables, manifests, or directories."""
    out: list[tuple[str | None, Path]] = []
    for p in paths:
        if p.is_dir():
            out.extend((None, f) for f in sorted(p.glob("*.tsv")) + sorted(p.glob("*/*.tsv")))
            continue
        if not p.is_file():
            print(f"WARNING: input not found: {p}", file=sys.stderr)
            continue
        first = ""
        with open(p, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                if raw.strip() and not raw.startswith("#"):
                    first = raw.strip()
                    break
        last_field = first.split("\t")[-1].strip() if first else ""
        if last_field and Path(last_field).is_file():
            out.extend((s, f) for s, f in read_manifest(p))
        else:
            out.append((None, p))
    return out


def parse_quast_report(path: Path) -> dict[str, dict[str, str]]:
    """Parse a QUAST report.tsv (metric rows x assembly columns)."""
    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines or lines[0].split("\t")[0].strip().lower() != "assembly":
        print(f"WARNING: {path} is not a QUAST report.tsv - skipping", file=sys.stderr)
        return out
    asm_cols = [clean_sample_name(c) for c in lines[0].split("\t")[1:]]
    for asm in asm_cols:
        out[asm] = {}
    for ln in lines[1:]:
        fields = ln.split("\t")
        metric = fields[0].strip()
        for asm, val in zip(asm_cols, fields[1:]):
            out[asm][metric] = val.strip()
    return out


def parse_checkm_table(path: Path) -> dict[str, dict[str, float | None]]:
    """Parse a CheckM qa / lineage_wf table."""
    out: dict[str, dict[str, float | None]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    hdr_idx = next((i for i, ln in enumerate(lines) if "Bin Id" in ln or "Bin ID" in ln), -1)
    if hdr_idx < 0:
        print(f"WARNING: no 'Bin Id' header in CheckM file {path} - skipping", file=sys.stderr)
        return out
    splitter = (lambda s: s.split("\t")) if "\t" in lines[hdr_idx] else (lambda s: re.split(r"\s{2,}", s.strip()))
    header = [h.strip() for h in splitter(lines[hdr_idx])]
    idx = {name: header.index(name) for name in ("Completeness", "Contamination", "Strain heterogeneity") if name in header}
    for ln in lines[hdr_idx + 1:]:
        if ln.lstrip().startswith(("-", "=")):
            continue
        parts = [p.strip() for p in splitter(ln)]
        if len(parts) < 3:
            continue

        def val(col: str) -> float | None:
            i = idx.get(col)
            return to_float(parts[i]) if i is not None and i < len(parts) else None

        out[clean_sample_name(parts[0])] = {
            "completeness": val("Completeness"),
            "contamination": val("Contamination"),
            "strain_het": val("Strain heterogeneity"),
        }
    return out


def parse_panaroo_matrix(
    path: Path,
    core: float,
    soft: float,
    shell: float,
    cohort: list[str],
) -> tuple[dict[str, dict[str, float]], int]:
    """Count per-isolate genes in each pangenome compartment from Panaroo output."""
    stats: dict[str, dict[str, float]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        delim = "\t" if "\t" in fh.readline() else ","
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        cols = [c for c in (reader.fieldnames or []) if c not in PANAROO_META_COLS]
        col_to_sample = {c: clean_sample_name(c) for c in cols}
        if cohort:
            keep = set(cohort)
            col_to_sample = {c: s for c, s in col_to_sample.items() if s in keep}
        n = len(col_to_sample)
        if n == 0:
            print(f"WARNING: no usable sample columns in Panaroo matrix {path}", file=sys.stderr)
            return stats, 0
        for s in col_to_sample.values():
            stats[s] = {"total": 0, "core": 0, "soft_core": 0, "shell": 0, "cloud": 0, "core_pct": 0.0}
        key_for = {"Core": "core", "SoftCore": "soft_core", "Shell": "shell", "Cloud": "cloud"}
        for row in reader:
            present = [s for c, s in col_to_sample.items()
                       if (row.get(c, "") or "").strip().lower() not in MISSING_TOKENS]
            if not present:
                continue
            comp = key_for[classify_compartment(len(present) / n, core, soft, shell)]
            for s in present:
                stats[s]["total"] += 1
                stats[s][comp] += 1
    for st in stats.values():
        st["core_pct"] = round(st["core"] / st["total"] * 100.0, 2) if st["total"] else 0.0
    return stats, n


def parse_mlst_table(path: Path) -> list[dict]:
    """Parse Torsten Seemann 'mlst' tabular output."""
    rows: list[dict] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3 or parts[0].upper() in ("FILE", "SAMPLE", "STRAIN"):
                continue
            rows.append({
                "sample": clean_sample_name(parts[0]),
                "scheme": parts[1].strip(),
                "st": format_st(parts[2]),
                "alleles": ";".join(p.strip() for p in parts[3:] if p.strip()),
            })
    return rows


def parse_ectyper_table(path: Path) -> list[dict]:
    """Parse ECTyper output.tsv."""
    rows: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for row in reader:
            name = pick(row, ("Name", "Sample", "sample"))
            if not name:
                continue
            o_type = pick(row, ("O-type", "O_type")) or "-"
            h_type = pick(row, ("H-type", "H_type")) or "-"
            sero = pick(row, ("Serotype", "serotype"))
            if not sero:
                sero = f"{o_type}:{h_type}" if o_type != "-" and h_type != "-" else o_type
            rows.append({
                "sample": clean_sample_name(name),
                "serotype": sero,
                "o_type": o_type,
                "h_type": h_type,
                "qc": pick(row, ("QC", "qc")),
            })
    return rows


def parse_kleborate_table(path: Path) -> list[dict]:
    """Parse Kleborate v2/v3 results TSV."""
    rows: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for row in reader:
            name = pick(row, ("strain", "Strain", "Sample", "sample"))
            if not name:
                continue
            rows.append({
                "sample": clean_sample_name(name),
                "st": format_st(pick(row, ("ST",), ("__mlst__ST", "__ST"))),
                "k_locus": pick(row, ("K_locus",), ("__K_locus",)) or "-",
                "k_type": pick(row, ("K_type",), ("__K_type",)) or "-",
                "o_locus": pick(row, ("O_locus",), ("__O_locus",)) or "-",
                "o_type": pick(row, ("O_type",), ("__O_type",)) or "-",
                "virulence_score": pick(row, ("virulence_score",), ("__virulence_score",)),
                "resistance_score": pick(row, ("resistance_score",), ("__resistance_score",)),
            })
    return rows


def parse_isescan(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse one isolate's ISEScan <genome>.fasta.csv."""
    if not base.exists():
        return None, []
    csv_path = first_file(base, ("*.fasta.csv", "*.csv"))
    if csv_path is None:
        return {"total": 0, "complete": 0, "partial": 0, "families": 0, "top_family": "-"}, []
    elements: list[dict] = []
    with open(csv_path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            t = (row.get("type", "") or "").strip().lower()
            status = "complete" if t in ("c", "complete") else "partial"
            elements.append({
                "Sample": sample,
                "Contig": row.get("seqID", "") or "",
                "IS_Family": row.get("family", "") or "Unknown",
                "IS_Cluster": row.get("cluster", "") or "",
                "Start": row.get("isBegin", "") or "",
                "End": row.get("isEnd", "") or "",
                "Length_bp": row.get("isLen", "") or "",
                "Strand": row.get("strand", "") or "",
                "Status": status,
                "Evalue": row.get("E-value", "") or "",
            })
    fams = Counter(e["IS_Family"] for e in elements)
    top = f"{fams.most_common(1)[0][0]} ({fams.most_common(1)[0][1]})" if fams else "-"
    summary = {
        "total": len(elements),
        "complete": sum(e["Status"] == "complete" for e in elements),
        "partial": sum(e["Status"] == "partial" for e in elements),
        "families": len(fams),
        "top_family": top,
    }
    return summary, elements


def parse_integron_finder(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse one isolate's IntegronFinder 2 output."""
    if not base.exists():
        return None, []
    integrons_file = first_file(base, ("*.integrons",))
    summary_file = first_file(base, ("*.summary",))
    empty = {"total": 0, "complete": 0, "in0": 0, "calin": 0, "attc": 0}

    if integrons_file is not None:
        groups: dict[tuple[str, str], dict] = {}
        with open(integrons_file, encoding="utf-8", errors="replace") as fh:
            lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
        if not lines or "ID_integron" not in lines[0]:
            return dict(empty), []
        for row in csv.DictReader(lines, delimiter="\t"):
            key = (row.get("ID_replicon", "") or "", row.get("ID_integron", "") or "")
            g = groups.setdefault(key, {
                "Sample": sample, "Integron_ID": key[1], "Contig": key[0],
                "Type": row.get("type", "") or "", "Start": None, "End": None,
                "attC_sites": 0, "Proteins": 0, "Has_intI": "no",
            })
            beg, end = to_int(row.get("pos_beg")), to_int(row.get("pos_end"))
            if beg is not None:
                g["Start"] = beg if g["Start"] is None else min(g["Start"], beg)
            if end is not None:
                g["End"] = end if g["End"] is None else max(g["End"], end)
            elt = (row.get("type_elt", "") or "").strip()
            if elt == "attC":
                g["attC_sites"] += 1
            elif elt == "protein":
                g["Proteins"] += 1
            elif elt == "intI":
                g["Has_intI"] = "yes"
            elif elt == "attI":
                g["Has_attI"] = "yes"
        records = list(groups.values())
        for r in records:
            r["Length_bp"] = (r["End"] - r["Start"] + 1) if r["Start"] is not None and r["End"] is not None else ""
        summary = {
            "total": len(records),
            "complete": sum(r["Type"].lower() == "complete" for r in records),
            "in0": sum(r["Type"].lower() == "in0" for r in records),
            "calin": sum(r["Type"].lower() == "calin" for r in records),
            "attc": sum(r["attC_sites"] for r in records),
        }
        return summary, records

    if summary_file is not None:
        totals = dict(empty)
        with open(summary_file, encoding="utf-8", errors="replace") as fh:
            lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
        for row in csv.DictReader(lines, delimiter="\t"):
            totals["complete"] += to_int(row.get("complete")) or 0
            totals["in0"] += to_int(row.get("In0")) or 0
            totals["calin"] += to_int(row.get("CALIN")) or 0
        totals["total"] = totals["complete"] + totals["in0"] + totals["calin"]
        totals["attc"] = ""
        return totals, []

    return dict(empty), []


def parse_phispy(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse one isolate's PhiSpy prophage coordinates."""
    if not base.exists():
        return None, []
    tsv = first_file(base, ("prophage_coordinates.tsv", "*prophage_coordinates.tsv", "prophage.tsv", "*prophage.tsv"))
    if tsv is None:
        return {"count": 0, "total_bp": 0, "longest_bp": 0}, []
    records: list[dict] = []
    with open(tsv, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if raw.startswith("#") or len(parts) < 4:
                continue
            start, stop = to_int(parts[2]), to_int(parts[3])
            if start is None or stop is None:
                continue
            records.append({
                "Sample": sample, "Prophage_ID": parts[0], "Contig": parts[1],
                "Start": start, "Stop": stop, "Length_bp": abs(stop - start) + 1,
            })
    lengths = [r["Length_bp"] for r in records]
    return {"count": len(records), "total_bp": sum(lengths), "longest_bp": max(lengths, default=0)}, records


def parse_platon(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """
    Parse Platon plasmid prediction outputs (<sample>.tsv).
    Header: ID, Length, Coverage, # RDS, Circular, Type ('plasmid' / 'chromosome').
    """
    if not base.exists():
        return None, []
    tsv_path = first_file(base, ("*.tsv",))
    if tsv_path is None:
        return {"plasmid_contigs": 0, "chromosome_contigs": 0, "plasmid_total_bp": 0, "circular_plasmids": 0}, []

    rows: list[dict] = []
    with open(tsv_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for r in reader:
            cid = pick(r, ("ID", "Contig", "contig_id", "seqID"))
            if not cid:
                continue
            ctype = pick(r, ("Type", "molecule_type", "type")).lower()
            clen = to_int(pick(r, ("Length", "length", "Length_bp", "len"))) or 0
            cov = pick(r, ("Coverage", "coverage", "depth"))
            circ = pick(r, ("Circular", "circular", "is_circular")).lower()
            is_circ = circ in ("yes", "true", "1")
            rows.append({
                "Sample": sample,
                "Contig": cid,
                "Length_bp": clen,
                "Coverage": cov,
                "Circular": "yes" if is_circ else "no",
                "Molecule_Type": "plasmid" if "plasmid" in ctype else "chromosome",
                "Source": "Platon",
            })

    p_rows = [r for r in rows if r["Molecule_Type"] == "plasmid"]
    c_rows = [r for r in rows if r["Molecule_Type"] == "chromosome"]
    summary = {
        "plasmid_contigs": len(p_rows),
        "chromosome_contigs": len(c_rows),
        "plasmid_total_bp": sum(r["Length_bp"] for r in p_rows),
        "circular_plasmids": sum(r["Circular"] == "yes" for r in p_rows),
    }
    return summary, rows


def parse_plasmidfinder(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse PlasmidFinder results (native results_tab.tsv or ABRicate plasmidfinder)."""
    if not base.exists():
        return None, []
    tsv_path = first_file(base, ("results_tab.tsv", "*results_tab.tsv", "*.tsv", "*.tab", "*.csv"))
    if tsv_path is None:
        return {"replicons_count": 0, "replicons": "-", "inc_groups": "-"}, []

    hits: list[dict] = []
    with open(tsv_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for r in reader:
            rep = pick(r, ("Plasmid", "GENE", "gene", "Replicon", "replicon"))
            if not rep:
                continue
            ident = pick(r, ("Identity", "%IDENTITY", "identity", "%Identity"))
            cov = pick(r, ("Query / Template length", "%COVERAGE", "coverage", "Coverage"))
            contig = pick(r, ("Contig", "SEQUENCE", "sequence", "contig_id"))
            acc = pick(r, ("Accession", "ACCESSION", "accession"))

            inc = "-"
            if rep.lower().startswith("inc"):
                m = re.match(r"^(inc[a-z0-9]+)", rep, re.IGNORECASE)
                inc = m.group(1).upper() if m else rep
            elif rep.lower().startswith("col"):
                inc = "Col"
            else:
                inc = rep.split("(")[0]

            hits.append({
                "Sample": sample,
                "Replicon": rep,
                "Inc_Group": inc,
                "Identity_%": ident,
                "Coverage_%": cov,
                "Contig": contig,
                "Accession": acc,
                "Source": "PlasmidFinder",
            })

    unique_reps = list(dict.fromkeys(h["Replicon"] for h in hits))
    unique_incs = list(dict.fromkeys(h["Inc_Group"] for h in hits if h["Inc_Group"] != "-"))
    unique_contigs = list(dict.fromkeys(h["Contig"] for h in hits if h["Contig"] and h["Contig"] != "-"))
    summary = {
        "replicons_count": len(hits),
        "replicons": ";".join(unique_reps) if unique_reps else "-",
        "inc_groups": ";".join(unique_incs) if unique_incs else "-",
        "plasmid_contigs": len(unique_contigs),
    }
    return summary, hits


def parse_mobsuite(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse MOB-suite outputs (contig_report.txt)."""
    if not base.exists():
        return None, []
    creport = first_file(base, ("contig_report.txt", "*contig_report.txt", "*contigs*.txt"))
    empty = {"plasmid_contigs": 0, "chromosome_contigs": 0, "plasmid_total_bp": 0,
             "replicons": "-", "predicted_mobility": "-"}
    if creport is None:
        return dict(empty), []

    records: list[dict] = []
    with open(creport, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for r in reader:
            cid = pick(r, ("contig_id", "Contig", "ID", "sequence_id"))
            if not cid:
                continue
            mtype = pick(r, ("molecule_type", "type")).lower()
            clen = to_int(pick(r, ("length", "size", "Length"))) or 0
            rep = pick(r, ("rep_type(s)", "rep_type", "rep_types"))
            relax = pick(r, ("relaxase_type(s)", "relaxase_type"))
            mpf = pick(r, ("mpf_type", "mpf_type(s)"))
            mob = pick(r, ("predicted_mobility", "mobility"))
            records.append({
                "Sample": sample,
                "Contig": cid,
                "Length_bp": clen,
                "Molecule_Type": "plasmid" if "plasmid" in mtype else "chromosome",
                "Replicon": rep if rep and rep not in MISSING_TOKENS else "-",
                "Relaxase": relax if relax and relax not in MISSING_TOKENS else "-",
                "MPF": mpf if mpf and mpf not in MISSING_TOKENS else "-",
                "Predicted_Mobility": mob if mob and mob not in MISSING_TOKENS else "-",
                "Source": "MOB-suite",
            })

    p_recs = [r for r in records if r["Molecule_Type"] == "plasmid"]
    c_recs = [r for r in records if r["Molecule_Type"] == "chromosome"]
    all_reps = [r["Replicon"] for r in p_recs if r["Replicon"] != "-"]
    all_mobs = [r["Predicted_Mobility"] for r in p_recs if r["Predicted_Mobility"] != "-"]
    summary = {
        "plasmid_contigs": len(p_recs),
        "chromosome_contigs": len(c_recs),
        "plasmid_total_bp": sum(r["Length_bp"] for r in p_recs),
        "replicons": ";".join(dict.fromkeys(all_reps)) if all_reps else "-",
        "predicted_mobility": ";".join(dict.fromkeys(all_mobs)) if all_mobs else "-",
    }
    return summary, records


def parse_breakpoints(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse detect_coverage_breakpoints.py outputs (breakpoints_report.tsv / .csv)."""
    if not base.exists():
        return None, []
    bfile = first_file(base, ("*breakpoints*.tsv", "*breakpoints*.csv", "breakpoints_report.tsv"))
    if bfile is None:
        return {"total": 0, "drops": 0, "jumps": 0, "gaps": 0}, []

    delim = "," if bfile.suffix.lower() == ".csv" else "\t"
    records: list[dict] = []
    with open(bfile, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter=delim)
        for r in reader:
            contig = pick(r, ("contig", "Contig", "seqID"))
            if not contig:
                continue
            btype = pick(r, ("type", "Type"))
            records.append({
                "Sample": sample,
                "Contig": contig,
                "Start": to_int(pick(r, ("start", "Start"))),
                "End": to_int(pick(r, ("end", "End"))),
                "Length_bp": to_int(pick(r, ("length_bp", "Length_bp", "length"))),
                "Type": btype,
                "Mean_Depth": pick(r, ("mean_depth", "Mean_Depth")),
                "Baseline_Depth": pick(r, ("baseline_depth", "Baseline_Depth")),
                "Baseline_Type": pick(r, ("baseline_type", "Baseline_Type")),
                "Fraction_or_Fold": pick(r, ("fraction_or_fold", "Fraction_or_Fold")),
                "Note": pick(r, ("note", "Note")),
            })

    summary = {
        "total": len(records),
        "drops": sum(r["Type"] == "coverage_drop" for r in records),
        "jumps": sum(r["Type"] == "coverage_jump_plateau" for r in records),
        "gaps": sum(r["Type"] == "zero_coverage_gap" for r in records),
    }
    return summary, records


def parse_annotation_summary(sample: str, base: Path) -> tuple[dict | None, list[dict]]:
    """Parse Prokka / Bakta summary file (<sample>.txt) or GFF to count rRNA, tRNA, repeat regions."""
    if not base.exists():
        return None, []

    txt_file = first_file(base, ("*.prokka.txt", "*_prokka.txt", "*.txt"))
    if txt_file is not None:
        rrna, trna, rep = 0, 0, 0
        with open(txt_file, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or ":" not in line:
                    continue
                k, _, v = line.partition(":")
                k = k.strip().lower()
                val = to_int(v.strip()) or 0
                if k in ("rrna", "rrnas"):
                    rrna += val
                elif k in ("trna", "trnas"):
                    trna += val
                elif k in ("repeat_region", "repeats", "repeat_regions"):
                    rep += val
        return {"rrna": rrna, "trna": trna, "repeat_region": rep}, []

    gff_file = first_file(base, ("*.gff3", "*.gff"))
    if gff_file is not None:
        rrna, trna, rep = 0, 0, 0
        with open(gff_file, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                if raw.startswith("#"):
                    continue
                parts = raw.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    ftype = parts[2].lower()
                    if ftype == "rrna":
                        rrna += 1
                    elif ftype == "trna":
                        trna += 1
                    elif ftype in ("repeat_region", "direct_repeat", "inverted_repeat"):
                        rep += 1
        return {"rrna": rrna, "trna": trna, "repeat_region": rep}, []

    return {"rrna": 0, "trna": 0, "repeat_region": 0}, []


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def load_typing(inputs: list[Path], parser) -> dict[str, dict]:
    """Load typing results (MLST / ECTyper / Kleborate) from tables or manifests."""
    out: dict[str, dict] = {}
    for manifest_sample, table in expand_table_inputs(inputs):
        if not table.is_file() or table.stat().st_size == 0:
            if manifest_sample:
                print(f"  [{manifest_sample}] WARNING: typing file missing/empty: {table}", file=sys.stderr)
            continue
        rows = parser(table)
        if manifest_sample and len(rows) == 1:
            rows[0]["sample"] = manifest_sample
        for r in rows:
            out[r["sample"]] = r
    return out


def load_quast(manifest: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    """Load QUAST metrics for every manifest entry."""
    data: dict[str, dict[str, str]] = {}
    order: list[str] = []
    for sample, path in read_manifest(manifest):
        report = first_file(path, ("report.tsv",))
        if report is None:
            print(f"  [{sample}] WARNING: QUAST report.tsv not found: {path}", file=sys.stderr)
            continue
        parsed = parse_quast_report(report)
        if len(parsed) == 1:
            parsed = {sample: next(iter(parsed.values()))}
        for s, metrics in parsed.items():
            data[s] = metrics
            for m in metrics:
                if m not in order:
                    order.append(m)
    return data, order


def load_manifest_data(manifest: Path, parser) -> tuple[dict[str, dict], list[dict]]:
    """Run a manifest parser over every entry ('sample<TAB>path')."""
    summaries: dict[str, dict] = {}
    elements: list[dict] = []
    for sample, path in read_manifest(manifest):
        summ, rows = parser(sample, path)
        if summ is None:
            print(f"  [{sample}] WARNING: output not found: {path}", file=sys.stderr)
            continue
        summaries[sample] = summ
        elements.extend(rows)
    return summaries, elements


def make_lookup(data: dict[str, dict]):
    """Return a function that finds a sample in data by exact, then normalised, name."""
    norm_index = {normalise(k): v for k, v in data.items()}

    def lookup(sample: str) -> dict:
        return data.get(sample) or norm_index.get(normalise(sample)) or {}

    return lookup


def evaluate_qc(
    rec: dict,
    species: str | None,
    min_comp: float | None,
    max_contam: float | None,
    max_contigs: int | None,
    min_core_pct: float | None,
    size_status: bool,
) -> dict[str, str]:
    """Compute the enabled QC status flags for one isolate."""
    flags: dict[str, str] = {}
    comp, contam = rec.get("_comp"), rec.get("_contam")
    if min_comp is not None or max_contam is not None:
        if comp is None or contam is None:
            flags["CheckM_Status"] = "UNKNOWN"
        elif comp >= (min_comp if min_comp is not None else 0.0) and contam <= (max_contam if max_contam is not None else 100.0):
            flags["CheckM_Status"] = "PASS"
        elif comp >= 90.0 and contam <= 10.0:
            flags["CheckM_Status"] = "WARNING"
        else:
            flags["CheckM_Status"] = "FAIL"

    if max_contigs is not None:
        n = rec.get("_contigs")
        flags["Assembly_Status"] = ("UNKNOWN" if n is None else "PASS" if n <= max_contigs
                                    else "WARNING" if n <= max_contigs + 200 else "FAIL")

    if size_status and species in SPECIES_SIZE_RANGES:
        lo, hi = SPECIES_SIZE_RANGES[species]
        L = rec.get("_length")
        flags["Size_Status"] = ("UNKNOWN" if L is None else "PASS" if lo <= L <= hi
                                else "WARNING" if lo * 0.9 <= L <= hi * 1.1 else "FAIL")

    if min_core_pct is not None:
        c = rec.get("_core_pct")
        flags["Pangenome_Status"] = ("UNKNOWN" if c is None else "PASS" if c >= min_core_pct
                                     else "WARNING" if c >= min_core_pct - 10.0 else "OUTLIER")

    if flags:
        known = [v for v in flags.values() if v != "UNKNOWN"]
        flags["Overall_QC"] = ("UNKNOWN" if not known else "FAIL" if any(v in ("FAIL", "OUTLIER") for v in known)
                               else "WARNING" if "WARNING" in known else "PASS")
    return flags


def build_overview(
    cohort: list[str],
    sources: dict[str, dict[str, dict]],
    species: str | None,
    qc_args: dict,
) -> tuple[list[dict], list[str], dict[str, list[str]]]:
    """Assemble one overview record per isolate."""
    seen: list[str] = []
    for data in sources.values():
        for s in data:
            if s not in seen:
                seen.append(s)
    if cohort:
        samples = list(cohort)
        cohort_norm = {normalise(s) for s in cohort}
        dropped = sorted(s for s in seen if normalise(s) not in cohort_norm)
    else:
        samples = sorted(seen)
        dropped = []

    lookups = {name: make_lookup(data) for name, data in sources.items()}
    missing: dict[str, list[str]] = {name: [] for name in sources}
    records: list[dict] = []

    for s in samples:
        got = {name: lk(s) for name, lk in lookups.items()}
        for name, d in got.items():
            if not d:
                missing[name].append(s)
        q, c, p = got.get("quast", {}), got.get("checkm", {}), got.get("panaroo", {})
        m, e, k = got.get("mlst", {}), got.get("ectyper", {}), got.get("kleborate", {})
        iss, itg, phi = got.get("isescan", {}), got.get("integron", {}), got.get("phispy", {})
        plt, mob = got.get("platon", {}), got.get("mobsuite", {})
        plf = got.get("plasmidfinder", {})
        bp = got.get("breakpoints", {})
        annot = got.get("annotation", {})

        length = to_int(q.get("Total length (>= 0 bp)") or q.get("Total length"))
        contigs = to_int(q.get("# contigs (>= 0 bp)") or q.get("# contigs"))
        gc = to_float(q.get("GC (%)"))

        # Chromosome vs Plasmid partitioning
        p_candidates = []
        if plt.get("plasmid_contigs") is not None and plt["plasmid_contigs"] > 0:
            p_candidates.append(plt["plasmid_contigs"])
        if mob.get("plasmid_contigs") is not None and mob["plasmid_contigs"] > 0:
            p_candidates.append(mob["plasmid_contigs"])
        if plf.get("plasmid_contigs") is not None and plf["plasmid_contigs"] > 0:
            p_candidates.append(plf["plasmid_contigs"])

        plasmid_contigs = max(p_candidates) if p_candidates else 0
        plasmid_bp = plt.get("plasmid_total_bp") or mob.get("plasmid_total_bp") or 0

        # Assembly total contigs from QUAST partitioned against plasmid contigs
        if contigs is not None:
            chr_contigs = max(0, contigs - plasmid_contigs)
        elif plt.get("chromosome_contigs") is not None and plt["chromosome_contigs"] > 0:
            chr_contigs = plt["chromosome_contigs"]
        elif mob.get("chromosome_contigs") is not None and mob["chromosome_contigs"] > 0:
            chr_contigs = mob["chromosome_contigs"]
        else:
            chr_contigs = None

        reps = plf.get("replicons") or mob.get("replicons") or ""
        inc = plf.get("inc_groups") or ""
        rrna = annot.get("rrna")
        trna = annot.get("trna")
        repeat_reg = annot.get("repeat_region")

        rec: dict = {
            "Sample": s,
            "Species": species or "",
            "_comp": c.get("completeness"),
            "_contam": c.get("contamination"),
            "_contigs": contigs,
            "_length": length,
            "_core_pct": p.get("core_pct"),
            "MLST_ST": m.get("st") or k.get("st") or ("-" if (m or k) else ""),
            "MLST_Scheme": m.get("scheme", ""),
            "MLST_Alleles": m.get("alleles", ""),
            "Serotype": e.get("serotype", ""),
            "Serotype_O": e.get("o_type", ""),
            "Serotype_H": e.get("h_type", ""),
            "ECTyper_QC": e.get("qc", ""),
            "K_Locus": k.get("k_locus", ""),
            "K_Type": k.get("k_type", ""),
            "O_Locus": k.get("o_locus", ""),
            "O_Type": k.get("o_type", ""),
            "Kleborate_Virulence_Score": k.get("virulence_score", ""),
            "Kleborate_Resistance_Score": k.get("resistance_score", ""),
            "MGE_IS_Total": iss.get("total", ""),
            "MGE_IS_Complete": iss.get("complete", ""),
            "MGE_Integron_Total": itg.get("total", ""),
            "MGE_Integron_Complete": itg.get("complete", ""),
            "MGE_Prophage_Count": phi.get("count", ""),
            "CheckM_Completeness_%": f"{c['completeness']:.2f}" if c.get("completeness") is not None else "",
            "CheckM_Contamination_%": f"{c['contamination']:.2f}" if c.get("contamination") is not None else "",
            "CheckM_Strain_Het": f"{c['strain_het']:.2f}" if c.get("strain_het") is not None else "",
            "QUAST_Total_Length_bp": length if length is not None else "",
            "QUAST_Contigs_Count": contigs if contigs is not None else "",
            "QUAST_N50_bp": q.get("N50", ""),
            "QUAST_GC_%": f"{gc:.2f}" if gc is not None else "",
            "Chromosome_Contigs_Count": chr_contigs if chr_contigs is not None else "",
            "Plasmid_Contigs_Count": plasmid_contigs if plasmid_contigs is not None else "",
            "Plasmid_Total_Length_bp": plasmid_bp if plasmid_bp is not None else "",
            "Plasmid_Replicons": reps if reps != "-" else "",
            "Plasmid_Inc_Groups": inc if inc != "-" else "",
            "rRNA_Count": rrna if rrna is not None else "",
            "tRNA_Count": trna if trna is not None else "",
            "Repeat_Region_Count": repeat_reg if repeat_reg is not None else "",
            "Breakpoint_Count": bp.get("total", ""),
            "Breakpoint_Drops": bp.get("drops", ""),
            "Breakpoint_Jumps": bp.get("jumps", ""),
            "Panaroo_Total_Genes": p.get("total", ""),
            "Panaroo_Core_Genes": p.get("core", ""),
            "Panaroo_Core_%": f"{p['core_pct']:.2f}" if p else "",
            "Panaroo_SoftCore_Genes": p.get("soft_core", ""),
            "Panaroo_Shell_Genes": p.get("shell", ""),
            "Panaroo_Cloud_Genes": p.get("cloud", ""),
        }
        rec.update(evaluate_qc(rec, species, **qc_args))
        records.append(rec)
    return records, dropped, missing


def select_columns(active: set[str], flags: list[str], species: str | None) -> list[str]:
    """Build the output column list from active sources and enabled QC flags."""
    cols = ["Sample"]
    if species:
        cols.append("Species")
    if flags:
        cols += ["Overall_QC"] + flags
    if "mlst" in active or "kleborate" in active:
        cols.append("MLST_ST")
    if "mlst" in active:
        cols += ["MLST_Scheme", "MLST_Alleles"]
    if "ectyper" in active:
        cols += ["Serotype", "Serotype_O", "Serotype_H", "ECTyper_QC"]
    if "kleborate" in active:
        cols += ["K_Locus", "K_Type", "O_Locus", "O_Type",
                 "Kleborate_Virulence_Score", "Kleborate_Resistance_Score"]
    if "isescan" in active:
        cols += ["MGE_IS_Total", "MGE_IS_Complete"]
    if "integron" in active:
        cols += ["MGE_Integron_Total", "MGE_Integron_Complete"]
    if "phispy" in active:
        cols += ["MGE_Prophage_Count"]
    if "checkm" in active:
        cols += ["CheckM_Completeness_%", "CheckM_Contamination_%", "CheckM_Strain_Het"]
    if "quast" in active:
        cols += ["QUAST_Total_Length_bp", "QUAST_Contigs_Count", "QUAST_N50_bp", "QUAST_GC_%"]
    if "platon" in active or "mobsuite" in active or "plasmidfinder" in active:
        cols += ["Chromosome_Contigs_Count", "Plasmid_Contigs_Count", "Plasmid_Total_Length_bp",
                 "Plasmid_Replicons", "Plasmid_Inc_Groups"]
    if "annotation" in active:
        cols += ["rRNA_Count", "tRNA_Count", "Repeat_Region_Count"]
    if "breakpoints" in active:
        cols += ["Breakpoint_Count", "Breakpoint_Drops", "Breakpoint_Jumps"]
    if "panaroo" in active:
        cols += ["Panaroo_Total_Genes", "Panaroo_Core_Genes", "Panaroo_Core_%",
                 "Panaroo_SoftCore_Genes", "Panaroo_Shell_Genes", "Panaroo_Cloud_Genes"]
    return cols


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

ISESCAN_COLUMNS = ["Sample", "Contig", "IS_Family", "IS_Cluster", "Start", "End",
                   "Length_bp", "Strand", "Status", "Evalue"]
INTEGRON_COLUMNS = ["Sample", "Integron_ID", "Contig", "Type", "Start", "End",
                    "Length_bp", "attC_sites", "Proteins", "Has_intI"]
PHISPY_COLUMNS = ["Sample", "Prophage_ID", "Contig", "Start", "Stop", "Length_bp"]
MGE_SUMMARY_COLUMNS = ["Sample", "IS_Total", "IS_Complete", "IS_Partial", "IS_Families",
                       "IS_Top_Family", "Integron_Total", "Integron_Complete", "Integron_In0",
                       "Integron_CALIN", "Integron_attC_sites", "Prophage_Count",
                       "Prophage_Total_bp", "Prophage_Longest_bp"]


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


def write_details(
    prefix: Path,
    cohort: list[str],
    quast: tuple[dict, list[str]] | None,
    mge: dict[str, tuple[dict, list[dict]]],
    plasmids: dict[str, tuple[dict, list[dict]]] | None = None,
    breakpoints: tuple[dict, list[dict]] | None = None,
) -> None:
    """Write full QUAST metrics, per-isolate MGE breakdown, plasmid tables and breakpoints."""
    keep = (lambda s: s in set(cohort)) if cohort else (lambda s: True)
    if quast:
        data, order = quast
        rows = [{"Sample": s, **data[s]} for s in sorted(data) if keep(s)]
        write_table(rows, ["Sample"] + order, Path(f"{prefix}_quast_full.tsv"),
                    ["TOOL: collate_genome_overview.py (all QUAST report.tsv metrics per isolate)"])
    if mge:
        samples = cohort or sorted({s for summ, _ in mge.values() for s in summ})
        iss = mge.get("isescan", ({}, []))[0]
        itg = mge.get("integron", ({}, []))[0]
        phi = mge.get("phispy", ({}, []))[0]
        rows = []
        for s in samples:
            a, b, c = iss.get(s, {}), itg.get(s, {}), phi.get(s, {})
            rows.append({
                "Sample": s,
                "IS_Total": a.get("total", ""), "IS_Complete": a.get("complete", ""),
                "IS_Partial": a.get("partial", ""), "IS_Families": a.get("families", ""),
                "IS_Top_Family": a.get("top_family", ""),
                "Integron_Total": b.get("total", ""), "Integron_Complete": b.get("complete", ""),
                "Integron_In0": b.get("in0", ""), "Integron_CALIN": b.get("calin", ""),
                "Integron_attC_sites": b.get("attc", ""),
                "Prophage_Count": c.get("count", ""), "Prophage_Total_bp": c.get("total_bp", ""),
                "Prophage_Longest_bp": c.get("longest_bp", ""),
            })
        write_table(rows, MGE_SUMMARY_COLUMNS, Path(f"{prefix}_mge_summary.tsv"), [
            "TOOL: collate_genome_overview.py (per-isolate MGE breakdown)",
            "IS: ISEScan type c=complete, p=partial. Integron: unique (contig, integron ID) from IntegronFinder .integrons; attC_sites approximates cassette count.",
            "Blank cell = tool output not found for that isolate (not the same as 0).",
        ])
        for name, cols, fname in (("isescan", ISESCAN_COLUMNS, "mge_isescan_elements"),
                                  ("integron", INTEGRON_COLUMNS, "mge_integrons"),
                                  ("phispy", PHISPY_COLUMNS, "mge_prophages")):
            if name in mge:
                rows = [r for r in mge[name][1] if keep(r["Sample"])]
                write_table(rows, cols, Path(f"{prefix}_{fname}.tsv"),
                            [f"TOOL: collate_genome_overview.py ({name} element-level table)"])

    if plasmids:
        pf_hits = []
        if "plasmidfinder" in plasmids:
            pf_hits.extend(plasmids["plasmidfinder"][1])
        if "mobsuite" in plasmids:
            for r in plasmids["mobsuite"][1]:
                if r.get("Replicon") and r["Replicon"] != "-":
                    pf_hits.append({
                        "Sample": r["Sample"],
                        "Replicon": r["Replicon"],
                        "Inc_Group": "-",
                        "Identity_%": "",
                        "Coverage_%": "",
                        "Contig": r["Contig"],
                        "Accession": "",
                        "Source": "MOB-suite",
                    })
        if pf_hits:
            reps_rows = [r for r in pf_hits if keep(r["Sample"])]
            write_table(reps_rows, ["Sample", "Replicon", "Inc_Group", "Identity_%", "Coverage_%", "Contig", "Accession", "Source"],
                        Path(f"{prefix}_plasmid_replicons.tsv"),
                        ["TOOL: collate_genome_overview.py (Plasmid replicon & Incompatibility group predictions)"])

        p_contigs = []
        if "platon" in plasmids:
            p_contigs.extend([r for r in plasmids["platon"][1] if r.get("Molecule_Type") == "plasmid"])
        if "mobsuite" in plasmids:
            p_contigs.extend([r for r in plasmids["mobsuite"][1] if r.get("Molecule_Type") == "plasmid"])
        if p_contigs:
            contig_rows = [r for r in p_contigs if keep(r["Sample"])]
            write_table(contig_rows, ["Sample", "Contig", "Length_bp", "Coverage", "Circular", "Molecule_Type", "Source"],
                        Path(f"{prefix}_plasmid_contigs.tsv"),
                        ["TOOL: collate_genome_overview.py (Predicted plasmid contigs)"])

    if breakpoints:
        bp_rows = [r for r in breakpoints[1] if keep(r["Sample"])]
        if bp_rows:
            write_table(bp_rows, ["Sample", "Contig", "Start", "End", "Length_bp", "Type", "Mean_Depth", "Baseline_Depth", "Baseline_Type", "Fraction_or_Fold", "Note"],
                        Path(f"{prefix}_coverage_breakpoints.tsv"),
                        ["TOOL: collate_genome_overview.py (Coverage breakpoints & depth anomalies)"])


def write_summary(
    records: list[dict],
    dropped: list[str],
    missing: dict[str, list[str]],
    flags: list[str],
    outpath: Path,
) -> None:
    """Write and echo the cohort audit summary."""
    n = len(records)
    pct = (lambda k: f"{k} ({k / n * 100:.1f}%)") if n else (lambda k: str(k))
    lines = [
        "=============================================================",
        " Master Genome Overview & QC Audit Summary",
        "=============================================================",
        f"  Isolates reported            : {n}",
        f"  Dropped (not in cohort list) : {len(dropped)} {dropped if dropped else ''}",
    ]
    for name, miss in missing.items():
        if miss:
            lines.append(f"  Missing {name:<12} data : {len(miss)} {miss}")
    lines.append("-------------------------------------------------------------")
    if flags:
        oq = Counter(r.get("Overall_QC", "") for r in records)
        for st in ("PASS", "WARNING", "FAIL", "UNKNOWN"):
            lines.append(f"  Overall_QC {st:<8}          : {pct(oq.get(st, 0))}")
        for f in flags:
            c = Counter(r.get(f, "") for r in records)
            lines.append(f"  {f:<17} PASS        : {pct(c.get('PASS', 0))}")
    else:
        lines.append("  QC flags                     : disabled (no thresholds given)")
    st = Counter(r["MLST_ST"] for r in records if r.get("MLST_ST") not in ("", "-"))
    sero = Counter(r["Serotype"] for r in records if r.get("Serotype") not in ("", "-"))
    lines += [
        "-------------------------------------------------------------",
        f"  Distinct STs                 : {len(st)}",
        f"  Top STs                      : {', '.join(f'{k} ({v})' for k, v in st.most_common(5)) or 'None'}",
        f"  Top serotypes                : {', '.join(f'{k} ({v})' for k, v in sero.most_common(5)) or 'None'}",
        "=============================================================",
    ]
    outpath.parent.mkdir(parents=True, exist_ok=True)
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
            "Collate QUAST, CheckM, Panaroo, MLST, ECTyper/Kleborate, MGEs, Plasmids, "
            "Coverage Breakpoints and Annotations into a master genome overview table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Required ---
    parser.add_argument(
        "--output", required=True, type=Path, metavar="TSV",
        help="Destination of the overview table; <output>.summary.txt is written alongside.",
    )

    # --- Cohort & inputs (all optional) ---
    parser.add_argument("--samples", type=Path, default=None, metavar="TXT",
                        help="Cohort list (one sample per line). Restricts and orders rows.")
    parser.add_argument("--quast-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample QUAST report.tsv files ('sample<TAB>path').")
    parser.add_argument("--checkm", type=Path, nargs="+", default=None, metavar="TSV",
                        help="One or more CheckM qa / lineage_wf tables.")
    parser.add_argument("--panaroo", type=Path, default=None, metavar="CSV",
                        help="Panaroo gene_presence_absence.csv (or .Rtab-like TSV).")
    parser.add_argument("--mlst", type=Path, nargs="+", default=None, metavar="TSV",
                        help="mlst tabular output(s) or manifests of per-sample mlst files.")
    parser.add_argument("--ectyper", type=Path, nargs="+", default=None, metavar="TSV",
                        help="ECTyper output.tsv file(s) or manifest (E. coli).")
    parser.add_argument("--kleborate", type=Path, nargs="+", default=None, metavar="TSV",
                        help="Kleborate v2/v3 results TSV(s) or manifest (Klebsiella).")
    parser.add_argument("--isescan-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample ISEScan output directories or .fasta.csv files.")
    parser.add_argument("--integron-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample IntegronFinder output directories (.integrons).")
    parser.add_argument("--phispy-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample PhiSpy output directories or prophage coordinate files.")

    # --- Plasmid, Breakpoint & Annotation manifests ---
    parser.add_argument("--platon-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample Platon output directories or TSV files ('sample<TAB>path').")
    parser.add_argument("--plasmidfinder-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample PlasmidFinder output tables ('sample<TAB>path').")
    parser.add_argument("--plasmidfinder", type=Path, nargs="+", default=None, metavar="TSV",
                        help="One or more PlasmidFinder / ABRicate plasmidfinder tables or manifests.")
    parser.add_argument("--mobsuite-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample MOB-suite directories or contig_report.txt files.")
    parser.add_argument("--breakpoints-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample detect_coverage_breakpoints.py outputs.")
    parser.add_argument("--annotation-manifest", type=Path, default=None, metavar="TXT",
                        help="Manifest of per-sample Prokka/Bakta summary txt or GFF files (counts rRNA, tRNA).")

    parser.add_argument("--detail-prefix", type=Path, default=None, metavar="PREFIX",
                        help="If set, also write <PREFIX>_quast_full.tsv, <PREFIX>_mge_summary.tsv, "
                             "_mge_*.tsv, _plasmid_*.tsv and _coverage_breakpoints.tsv.")

    # --- QC thresholds (flags disabled unless set) ---
    parser.add_argument("--species", choices=sorted(SPECIES_SIZE_RANGES), default=None,
                        help="Species key; adds a Species column and enables --size-status ranges.")
    parser.add_argument("--size-status", action="store_true",
                        help="Add Size_Status using the expected genome-size range of --species.")
    parser.add_argument("--min-completeness", type=float, default=None, metavar="FLOAT",
                        help="CheckM completeness %% required for PASS (e.g. 95). Enables CheckM_Status.")
    parser.add_argument("--max-contamination", type=float, default=None, metavar="FLOAT",
                        help="CheckM contamination %% allowed for PASS (e.g. 5). Enables CheckM_Status.")
    parser.add_argument("--max-contigs", type=int, default=None, metavar="INT",
                        help="Maximum contig count for PASS (e.g. 300). Enables Assembly_Status.")
    parser.add_argument("--min-core-pct", type=float, default=None, metavar="FLOAT",
                        help="Minimum %% of an isolate's genes that are Core for PASS (e.g. 80).")
    parser.add_argument("--core-cutoff", type=float, default=0.99, metavar="FLOAT",
                        help="Cohort frequency for Core genes. Default: 0.99.")
    parser.add_argument("--soft-core-cutoff", type=float, default=0.95, metavar="FLOAT",
                        help="Cohort frequency for Soft-core genes. Default: 0.95.")
    parser.add_argument("--shell-cutoff", type=float, default=0.15, metavar="FLOAT",
                        help="Cohort frequency for Shell genes (below = Cloud). Default: 0.15.")

    return parser.parse_args()


def main():
    """Orchestrate loading of every source, collation and writing."""
    args = parse_args()

    inputs_given = [
        args.quast_manifest, args.checkm, args.panaroo, args.mlst, args.ectyper,
        args.kleborate, args.isescan_manifest, args.integron_manifest, args.phispy_manifest,
        args.platon_manifest, args.plasmidfinder_manifest, args.plasmidfinder,
        args.mobsuite_manifest, args.breakpoints_manifest, args.annotation_manifest
    ]
    if not any(inputs_given):
        sys.exit("ERROR: provide at least one input source.")
    if args.size_status and not args.species:
        sys.exit("ERROR: --size-status requires --species.")

    cohort = read_sample_list(args.samples)
    qc_args = {
        "min_comp": args.min_completeness, "max_contam": args.max_contamination,
        "max_contigs": args.max_contigs, "min_core_pct": args.min_core_pct,
        "size_status": args.size_status,
    }

    print("\n========================================================")
    print(" collate_genome_overview.py - master genome overview")
    print("========================================================")
    print(f" Cohort list : {args.samples or 'None (union of all sources)'}"
          f"{f' ({len(cohort)} samples)' if cohort else ''}")
    for label, val in (
        ("QUAST", args.quast_manifest), ("CheckM", args.checkm), ("Panaroo", args.panaroo),
        ("MLST", args.mlst), ("ECTyper", args.ectyper), ("Kleborate", args.kleborate),
        ("ISEScan", args.isescan_manifest), ("Integron", args.integron_manifest),
        ("PhiSpy", args.phispy_manifest), ("Platon", args.platon_manifest),
        ("PlasmidFinder", args.plasmidfinder_manifest or args.plasmidfinder),
        ("MOB-suite", args.mobsuite_manifest), ("Breakpoints", args.breakpoints_manifest),
        ("Annotation", args.annotation_manifest)
    ):
        if val:
            print(f" {label:<14} : {val}")
    print(f" Output         : {args.output}")
    print("========================================================\n")

    sources: dict[str, dict[str, dict]] = {}
    quast_full = None
    mge: dict[str, tuple[dict, list[dict]]] = {}
    plasmids_detail: dict[str, tuple[dict, list[dict]]] = {}
    breakpoints_detail: tuple[dict, list[dict]] | None = None

    if args.quast_manifest:
        quast_full = load_quast(args.quast_manifest)
        sources["quast"] = quast_full[0]
    if args.checkm:
        checkm: dict[str, dict] = {}
        for t in args.checkm:
            if t.is_file():
                checkm.update(parse_checkm_table(t))
            else:
                print(f"WARNING: CheckM table not found: {t}", file=sys.stderr)
        sources["checkm"] = checkm
    if args.panaroo:
        if args.panaroo.is_file():
            pan, n_gen = parse_panaroo_matrix(args.panaroo, args.core_cutoff, args.soft_core_cutoff,
                                              args.shell_cutoff, cohort)
            print(f"Panaroo: compartments computed over {n_gen} genomes.")
            sources["panaroo"] = pan
        else:
            print(f"WARNING: Panaroo matrix not found: {args.panaroo}", file=sys.stderr)
            sources["panaroo"] = {}
    if args.mlst:
        sources["mlst"] = load_typing(args.mlst, parse_mlst_table)
    if args.ectyper:
        sources["ectyper"] = load_typing(args.ectyper, parse_ectyper_table)
    if args.kleborate:
        sources["kleborate"] = load_typing(args.kleborate, parse_kleborate_table)

    for name, manifest, fn in (("isescan", args.isescan_manifest, parse_isescan),
                               ("integron", args.integron_manifest, parse_integron_finder),
                               ("phispy", args.phispy_manifest, parse_phispy)):
        if manifest:
            mge[name] = load_manifest_data(manifest, fn)
            sources[name] = mge[name][0]

    # Plasmids, Breakpoints & Annotation
    if args.platon_manifest:
        plt_data = load_manifest_data(args.platon_manifest, parse_platon)
        sources["platon"] = plt_data[0]
        plasmids_detail["platon"] = plt_data

    if args.plasmidfinder_manifest:
        pf_data = load_manifest_data(args.plasmidfinder_manifest, parse_plasmidfinder)
        sources["plasmidfinder"] = pf_data[0]
        plasmids_detail["plasmidfinder"] = pf_data
    elif args.plasmidfinder:
        pf_loaded = load_typing(args.plasmidfinder, lambda p: parse_plasmidfinder(sample_from_path(p), p)[1])
        pf_summ: dict[str, dict] = {}
        all_pf_hits = []
        for s, r_list in pf_loaded.items():
            reps = list(dict.fromkeys(r.get("Replicon") for r in r_list if r.get("Replicon")))
            incs = list(dict.fromkeys(r.get("Inc_Group") for r in r_list if r.get("Inc_Group") and r.get("Inc_Group") != "-"))
            pf_summ[s] = {"replicons_count": len(reps), "replicons": ";".join(reps) if reps else "-", "inc_groups": ";".join(incs) if incs else "-"}
            all_pf_hits.extend(r_list)
        sources["plasmidfinder"] = pf_summ
        plasmids_detail["plasmidfinder"] = (pf_summ, all_pf_hits)

    if args.mobsuite_manifest:
        mob_data = load_manifest_data(args.mobsuite_manifest, parse_mobsuite)
        sources["mobsuite"] = mob_data[0]
        plasmids_detail["mobsuite"] = mob_data

    if args.breakpoints_manifest:
        bp_data = load_manifest_data(args.breakpoints_manifest, parse_breakpoints)
        sources["breakpoints"] = bp_data[0]
        breakpoints_detail = bp_data

    if args.annotation_manifest:
        annot_data = load_manifest_data(args.annotation_manifest, parse_annotation_summary)
        sources["annotation"] = annot_data[0]

    records, dropped, missing = build_overview(cohort, sources, args.species, qc_args)
    if not records:
        sys.exit("ERROR: no overview records produced - check inputs and --samples.")

    flags = [f for f, on in (
        ("CheckM_Status", args.min_completeness is not None or args.max_contamination is not None),
        ("Assembly_Status", args.max_contigs is not None),
        ("Size_Status", args.size_status),
        ("Pangenome_Status", args.min_core_pct is not None),
    ) if on]
    columns = select_columns(set(sources), flags, args.species)

    comments = [
        "TOOL: collate_genome_overview.py (Master Bacterial Genome Overview & QC)",
        f"SOURCES: {', '.join(sorted(sources))}",
        f"COHORT: {len(records)} isolates" + (f" from {args.samples.name}; {len(dropped)} non-cohort isolates dropped"
                                              if cohort else " (union of all sources)"),
    ]
    if flags:
        th = []
        if "CheckM_Status" in flags:
            th.append(f"--min-completeness {args.min_completeness} --max-contamination {args.max_contamination}")
        if "Assembly_Status" in flags:
            th.append(f"--max-contigs {args.max_contigs}")
        if "Size_Status" in flags:
            lo, hi = SPECIES_SIZE_RANGES[args.species]
            th.append(f"--species {args.species} --size-status")
        if "Pangenome_Status" in flags:
            th.append(f"--min-core-pct {args.min_core_pct}")
        comments.append(f"THRESHOLDS: {' '.join(th)}")
        if "CheckM_Status" in flags:
            comments.append(f"  CheckM_Status: PASS comp >= {args.min_completeness or 0.0}% and contam <= "
                            f"{args.max_contamination if args.max_contamination is not None else 100.0}%; "
                            "WARNING comp >= 90% and contam <= 10%; else FAIL")
        if "Assembly_Status" in flags:
            comments.append(f"  Assembly_Status: PASS contigs <= {args.max_contigs}; WARNING <= {args.max_contigs + 200}; else FAIL")
        if "Size_Status" in flags:
            comments.append(f"  Size_Status: expected {args.species} genome {lo / 1e6:.2f}-{hi / 1e6:.2f} Mb = PASS; within +/-10% = WARNING; else FAIL")
        if "Pangenome_Status" in flags:
            comments.append(f"  Pangenome_Status: PASS Core_% >= {args.min_core_pct}; WARNING >= {args.min_core_pct - 10.0}; else OUTLIER")
        comments.append("  Overall_QC: FAIL if any FAIL/OUTLIER; WARNING if any WARNING; else PASS (UNKNOWN = data missing)")
    else:
        comments.append("THRESHOLDS: none (QC flags disabled; raw metrics only)")

    if "checkm" in sources:
        comments.append("CheckM_Strain_Het: CheckM 'Strain heterogeneity' = % of multi-copy markers with >= 90% AAI between copies.")
    if "panaroo" in sources:
        comments.append(f"PANAROO: Core >= {args.core_cutoff * 100:.0f}%, SoftCore >= {args.soft_core_cutoff * 100:.0f}%, "
                        f"Shell >= {args.shell_cutoff * 100:.0f}%, Cloud < {args.shell_cutoff * 100:.0f}% of cohort genomes.")
    if mge:
        comments.append("MGE: IS = ISEScan (Complete = type c); Integron = unique IntegronFinder integrons "
                        "(Complete = attC + intI); Prophage = PhiSpy regions. Blank = tool output not found.")
    if plasmids_detail:
        comments.append("PLASMIDS: Platon / PlasmidFinder / MOB-suite plasmid and replicon predictions.")
    if breakpoints_detail:
        comments.append("BREAKPOINTS: detect_coverage_breakpoints.py depth drops & jumps.")
    if "annotation" in sources:
        comments.append("ANNOTATION: Prokka / Bakta rRNA, tRNA and repeat region counts.")

    write_table(records, columns, args.output, comments)
    write_summary(records, dropped, missing, flags, args.output.with_suffix(".summary.txt"))
    if args.detail_prefix:
        write_details(args.detail_prefix, cohort, quast_full, mge, plasmids_detail, breakpoints_detail)

    print("\nDone.")


if __name__ == "__main__":
    main()

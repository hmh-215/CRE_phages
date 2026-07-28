"""
ncbi_blastn_auto.py
════════════════════════════════════════════════════════════════════
Automates remote BLASTn searches against NCBI's full nr database
using Biopython's NCBIWWW / NCBIXML modules.

Each sequence in the input FASTA is submitted as a SEPARATE job.
The job name is automatically derived from the FASTA header as:
    <sample_name>_<node_number>
e.g., header ">NODE_5_length_3200_cov_45.2" → job name "NODE_5"

Output: one TSV file per sequence  (BLAST -outfmt 6 style columns)
        one summary TSV listing all job names and their output files

Usage
─────
    python ncbi_blastn_auto.py \
        --input   /path/to/sequences.fasta   \
        --output  /path/to/results/          \
        --email   your_email@example.com     \
        --hits    10

Optional flags
─────────────
    --hits   N        top N hits to return (default: 10)
    --apikey KEY      NCBI API key (raises rate limit 3→10 req/s)
    --delay  N        seconds to wait between submissions (default: 5)

Requirements
────────────
    pip install biopython
"""

import argparse
import os
import re
import time
from datetime import datetime
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML


# ──────────────────────────────────────────────────────────────────
# Job-name extraction
# ──────────────────────────────────────────────────────────────────

def derive_job_name(header: str) -> str:
    """
    Build a job name from a FASTA header.

    Strategy (tries each pattern in order):
      1. SPAdes-style: 'NODE_5_length_...' → 'NODE_5'
      2. Any 'NODE_<digits>' token          → 'NODE_<digits>'
      3. First two underscore-separated tokens from the seq ID
      4. Fallback: the raw seq ID, sanitised

    Returns a filesystem-safe string.
    """
    seq_id = header.split()[0].lstrip(">")

    # Pattern 1 – SPAdes header (most common in phage/metagenome assemblies)
    m = re.match(r"(NODE_\d+)", seq_id, re.IGNORECASE)
    if m:
        return m.group(1)

    # Pattern 2 – generic "NODE_N" anywhere in the header
    m = re.search(r"(NODE_\d+)", header, re.IGNORECASE)
    if m:
        return m.group(1)

    # Pattern 3 – take first two tokens split on underscore
    parts = seq_id.split("_")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"

    # Fallback – sanitise the raw ID
    safe = re.sub(r"[^\w\-.]", "_", seq_id)
    return safe[:60]   # cap length


# ──────────────────────────────────────────────────────────────────
# BLAST submission
# ──────────────────────────────────────────────────────────────────

def submit_blast(sequence: str, hits: int) -> object:
    """
    Submit one nucleotide sequence to NCBI BLASTn (nr database).
    Returns an open file-like handle containing the XML result.

    Key qblast parameters
    ─────────────────────
    program      'blastn'   – nucleotide vs nucleotide
    database     'nr'       – full non-redundant NT database
    hitlist_size              max hits per query
    format_type  'XML'      – required for Biopython's NCBIXML parser
    """
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nr",
        sequence=sequence,
        hitlist_size=hits,
        format_type="XML",
    )
    return result_handle


# ──────────────────────────────────────────────────────────────────
# Result parsing & saving
# ──────────────────────────────────────────────────────────────────

# Output columns match BLAST -outfmt 6
TSV_HEADER = (
    "job_name\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\t"
    "qstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\n"
)


def parse_and_save(result_handle, output_dir: str, job_name: str, hits: int) -> str:
    """
    Parse the BLAST XML result and write a tab-delimited file.
    Returns the path of the output file.
    """
    blast_records = NCBIXML.parse(result_handle)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_file = os.path.join(output_dir, f"{job_name}_{timestamp}.tsv")

    rows_written = 0
    with open(out_file, "w") as fh:
        fh.write(TSV_HEADER)
        for blast_record in blast_records:
            q_id = blast_record.query.split()[0]
            for alignment in blast_record.alignments[:hits]:
                for hsp in alignment.hsps[:1]:          # best HSP only
                    align_len = hsp.align_length or 1   # avoid /0
                    pident = round(hsp.identities / align_len * 100, 2)
                    mismatches = align_len - hsp.identities - hsp.gaps

                    row = (
                        f"{job_name}\t"
                        f"{q_id}\t"
                        f"{alignment.accession}\t"
                        f"{pident}\t"
                        f"{align_len}\t"
                        f"{mismatches}\t"
                        f"{hsp.gaps}\t"
                        f"{hsp.query_start}\t{hsp.query_end}\t"
                        f"{hsp.sbjct_start}\t{hsp.sbjct_end}\t"
                        f"{hsp.expect:.2e}\t"
                        f"{hsp.bits:.1f}\t"
                        f"{alignment.title[:150]}\n"
                    )
                    fh.write(row)
                    rows_written += 1

    print(f"        → {rows_written} hit(s) saved to: {os.path.basename(out_file)}")
    return out_file


# ──────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Remote BLASTn (nr) via Biopython — auto job-naming from FASTA headers."
    )
    parser.add_argument("--input",  required=True,
                        help="Path to input FASTA file (multi-sequence supported)")
    parser.add_argument("--output", required=True,
                        help="Directory where result TSV files will be saved")
    parser.add_argument("--email",  required=True,
                        help="Your e-mail address (mandatory per NCBI Entrez policy)")
    parser.add_argument("--hits",   type=int, default=10,
                        help="Number of top hits per query (default: 10)")
    parser.add_argument("--apikey", default=None,
                        help="NCBI API key (optional; raises rate limit from 3 to 10 req/s)")
    parser.add_argument("--delay",  type=float, default=5.0,
                        help="Seconds to wait between submissions (default: 5)")
    args = parser.parse_args()

    # ── Validate ──────────────────────────────────────────────────
    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    os.makedirs(args.output, exist_ok=True)

    # ── Configure Entrez ──────────────────────────────────────────
    Entrez.email = args.email
    if args.apikey:
        Entrez.api_key = args.apikey
        print(f"[INFO] NCBI API key loaded (rate limit: 10 req/s)")
    else:
        print(f"[INFO] No API key — rate limit: 3 req/s (delay {args.delay}s between jobs)")

    # ── Load sequences ────────────────────────────────────────────
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {args.input}. Check the file is valid FASTA.")
    print(f"[INFO] Loaded {len(records)} sequence(s) from {args.input}\n")

    # ── Summary log ───────────────────────────────────────────────
    summary_path = os.path.join(
        args.output,
        f"blast_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"
    )
    summary_lines = ["seq_index\tjob_name\tseq_id\tseq_length\tresult_file\tstatus\n"]

    # ── Main loop ─────────────────────────────────────────────────
    for i, record in enumerate(records, start=1):
        job_name = derive_job_name(record.description)
        seq_len  = len(record.seq)

        print(f"[{i:>3}/{len(records)}]  Job: {job_name}  |  ID: {record.id}  |  Length: {seq_len} bp")

        try:
            print(f"        Submitting to NCBI BLASTn ... (may take several minutes)")
            result_handle = submit_blast(str(record.seq), args.hits)
            out_file = parse_and_save(result_handle, args.output, job_name, args.hits)
            summary_lines.append(
                f"{i}\t{job_name}\t{record.id}\t{seq_len}\t{os.path.basename(out_file)}\tOK\n"
            )

        except Exception as exc:
            print(f"        [ERROR] {exc}")
            summary_lines.append(
                f"{i}\t{job_name}\t{record.id}\t{seq_len}\tN/A\tERROR: {exc}\n"
            )

        # Polite delay — NCBI blocks IPs that submit too fast
        if i < len(records):
            print(f"        Waiting {args.delay}s before next submission …\n")
            time.sleep(args.delay)

    # ── Write summary ─────────────────────────────────────────────
    with open(summary_path, "w") as sf:
        sf.writelines(summary_lines)

    print(f"\n{'─'*60}")
    print(f"[DONE]  {len(records)} sequences processed.")
    print(f"        Results directory : {args.output}")
    print(f"        Summary log       : {summary_path}")
    print(f"{'─'*60}")


if __name__ == "__main__":
    main()

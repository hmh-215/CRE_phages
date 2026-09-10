#!/bin/bash
# script for annotating bacteriophage ILLUMINA samples (Nithesis-like phage)
# Building No. 1
set -euo pipefail

################
# GLOBAL SETUP #
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/nithesis"
MAPPING_PHAGES="${SAMPLE_PATH}/phages/mapping_phages"
threads=16

# Paths to databases
checkvdb="${REF_PATH}/checkv-db-v1.5"
eggnog_db="${REF_PATH}/eggnog_db"
pharokka_db="${REF_PATH}/pharokka_db"

# Create working directories
STRUCTURE_NITHESIS="${WORK_PATH}/structure_nithesis"
ANNOTATION_NITHESIS="${WORK_PATH}/annotation_nithesis"
PHAGETERM_REF="${WORK_PATH}/phageterm_references"

mkdir -p "${STRUCTURE_NITHESIS}"
mkdir -p "${ANNOTATION_NITHESIS}"
mkdir -p "${PHAGETERM_REF}"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/nithesis_annotation.log"
FAIL_LOG="${WORK_PATH}/nithesis_annotation.failed_skipped.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================" >> "${FAIL_LOG}"
echo " Failure & Skip Log â€” Started: $(date)" >> "${FAIL_LOG}"
echo "======================================================" >> "${FAIL_LOG}"

# Helper function to record skipped or failed steps with timestamp
log_failure() {
    local phase="$1"
    local sample="$2"
    local status="$3" # e.g. "SKIPPED_EXISTS", "INPUT_MISSING", "EXECUTION_FAILED", "OUTPUT_MISSING_OR_EMPTY"
    local reason="$4"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${phase}] [${sample}] [${status}] ${reason}" >> "${FAIL_LOG}"
}

echo "============================================================"
echo " Nithesis-like phage annotations - Started: $(date)"
echo " Full log               : ${LOG}"
echo " Failure & skip log     : ${FAIL_LOG}"
echo "============================================================"

    #################################################
    # PHASE 1: CHECK GENOME COMPLETENESS (CHECKV)   #
    #################################################

for sample_id in 22 23 24; do
    sample_name="WS2762512A${sample_id}"
    checkv_out_dir="${STRUCTURE_NITHESIS}/${sample_name}.checkv"
    checkv_summary="${checkv_out_dir}/quality_summary.tsv"
    filtered_assembly="${SAMPLE_PATH}/phages/assembly_phages/${sample_name}.contigs.filtered.fasta"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m CHECKV: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${checkv_summary}" ]; then
        echo -e "\e[32m   [${sample_name}] CheckV output already exists -> skipping \e[0m"
        log_failure "Phase1_CheckV" "${sample_name}" "SKIPPED_EXISTS" "Output ${checkv_summary} already exists"
    elif [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Assembly missing: ${filtered_assembly} -> skipping CheckV \e[0m"
        log_failure "Phase1_CheckV" "${sample_name}" "INPUT_MISSING" "Assembly ${filtered_assembly} missing"
    else
        rm -rf "${checkv_out_dir}"
        if ! conda run -n BPstructure checkv end_to_end \
            -t "${threads}" \
            -d "${checkvdb}" \
            "${filtered_assembly}" \
            "${checkv_out_dir}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: CheckV failed \e[0m"
            log_failure "Phase1_CheckV" "${sample_name}" "EXECUTION_FAILED" "checkv non-zero exit status"
        else
            echo -e "\e[32m CheckV complete for ${sample_name} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 2: EXTRACT HIGH-COMPLETENESS CONTIGS    #
    #################################################

    echo -e "\e[32m Extracting high-completeness contigs (>=90%) from CheckV results... \e[0m"

for sample_id in 21 22 23 24; do
    sample_name="WS2762512A${sample_id}"
    checkv_tsv="${STRUCTURE_NITHESIS}/${sample_name}.checkv/quality_summary.tsv"
    filtered_assembly="${SAMPLE_PATH}/phages/assembly_phages/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${checkv_tsv}" ]; then
        echo -e "\e[33m CheckV summary not found for ${sample_name} (${checkv_tsv}), skipping node extraction \e[0m"
        continue
    fi
    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[33m Filtered assembly not found for ${sample_name}, skipping node extraction \e[0m"
        continue
    fi

    # Parse TSV: col 1 = contig_id, col 9 or 10 = completeness (%)
    tail -n +2 "${checkv_tsv}" | awk -F'\t' '($9 != "NA" && $9+0 >= 90) || ($10 != "NA" && $10+0 >= 90) { print $1 }' | \
    while read -r contig_id; do
        node_num=$(echo "${contig_id}" | grep -oP '(?<=NODE_)\d+' || true)
        if [ -z "${node_num}" ]; then
            continue
        fi

        out_fasta="${PHAGETERM_REF}/${sample_name}_NODE_${node_num}.metaSPAdes.fasta"
        new_header="${sample_name}_NODE_${node_num}.metaSPAdes"

        if [ -s "${out_fasta}" ]; then
            echo -e "\e[32m   Node ${node_num} already extracted for ${sample_name} -> skipping \e[0m"
            log_failure "Phase2_NodeExtract" "${sample_name}" "SKIPPED_EXISTS" "File ${out_fasta} already exists"
        else
            conda run -n assembly seqkit grep \
                -p "${contig_id}" \
                "${filtered_assembly}" \
                | conda run -n assembly seqkit replace \
                --pattern ".*" \
                --replacement "${new_header}" \
                > "${out_fasta}"
            echo -e "\e[32m   Extracted: ${out_fasta} (contig: ${contig_id}) \e[0m"
        fi
    done
done

    #################################################
    # PHASE 3: PREDICTION OF PHAGE TERMINI          #
    #################################################

for phageterm_ref in "${PHAGETERM_REF}"/*.metaSPAdes.fasta; do
    [ -f "${phageterm_ref}" ] || continue

    basename_no_ext=$(basename "${phageterm_ref}" .metaSPAdes.fasta)
    sample_name=$(echo "${basename_no_ext}" | grep -oP '^WS\d+A\d+' || true)
    node_num=$(echo "${basename_no_ext}" | grep -oP '(?<=NODE_)\d+' || true)

    if [ -z "${sample_name}" ] || [ -z "${node_num}" ]; then
        continue
    fi

    TERM_OUT="${STRUCTURE_NITHESIS}/${sample_name}_NODE_${node_num}.term"
    mkdir -p "${TERM_OUT}"

    read1f="${MAPPING_PHAGES}/${sample_name}/${sample_name}.R1.f12.fastq.gz"
    read2f="${MAPPING_PHAGES}/${sample_name}/${sample_name}.R2.f12.fastq.gz"

    if [ ! -s "${read1f}" ] || [ ! -s "${read2f}" ]; then
        # Fallback to WORK_PATH if mapped in nithesis folder
        read1f="${WORK_PATH}/mapping_phages/${sample_name}/${sample_name}.R1.f12.fastq.gz"
        read2f="${WORK_PATH}/mapping_phages/${sample_name}/${sample_name}.R2.f12.fastq.gz"
    fi

    if [ ! -s "${read1f}" ] || [ ! -s "${read2f}" ]; then
        echo -e "\e[31m   [${sample_name}_NODE_${node_num}] Filtered reads missing -> skipping PhageTerm \e[0m"
        log_failure "Phase3_PhageTerm" "${sample_name}_NODE_${node_num}" "INPUT_MISSING" "Filtered reads missing"
        continue
    fi

    echo -e "\e[31m ======================== \e[0m"
    echo -e "\e[31m PHAGETERM: ${sample_name} NODE_${node_num} \e[0m"
    echo -e "\e[31m ======================== \e[0m"

    expected_term_report="${TERM_OUT}/${sample_name}_NODE_${node_num}_report.pdf"
    if [ -s "${expected_term_report}" ]; then
        echo -e "\e[32m   [${sample_name}_NODE_${node_num}] PhageTerm report already exists -> skipping \e[0m"
        log_failure "Phase3_PhageTerm" "${sample_name}_NODE_${node_num}" "SKIPPED_EXISTS" "Output ${expected_term_report} already exists"
    else
        if ! conda run -n phageterm_env phageterm \
            -f "${read1f}" \
            -p "${read2f}" \
            -r "${phageterm_ref}" \
            -c "${threads}" \
            -o "${TERM_OUT}" \
            --report_title "${sample_name}_NODE_${node_num}"; then
            echo -e "\e[31m   [${sample_name}_NODE_${node_num}] ERROR: PhageTerm failed \e[0m"
            log_failure "Phase3_PhageTerm" "${sample_name}_NODE_${node_num}" "EXECUTION_FAILED" "phageterm non-zero exit status"
        else
            echo -e "\e[32m PhageTerm complete for ${sample_name}_NODE_${node_num} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 4: MULTIPLE SEQUENCE ALIGNMENT (MAFFT)  #
    #################################################

    mafft_input="${ANNOTATION_NITHESIS}/mafft_input_nithesis.fasta"
    mafft_output="${ANNOTATION_NITHESIS}/mafft_output_nithesis.aln"

    if [ -s "${mafft_output}" ]; then
        echo -e "\e[32m MAFFT alignment output already exists -> skipping \e[0m"
        log_failure "Phase4_MAFFT" "ALL" "SKIPPED_EXISTS" "Output ${mafft_output} already exists"
    elif [ -s "${mafft_input}" ]; then
        echo -e "\e[31m ============================= \e[0m"
        echo -e "\e[31m MAFFT: ALIGN NITHESIS CONTIGS \e[0m"
        echo -e "\e[31m ============================= \e[0m"

        conda run -n BPannotate mafft \
            --thread "${threads}" \
            --clustalout \
            "${mafft_input}" > "${mafft_output}" 2>/dev/null || true
        echo -e "\e[32m MAFFT alignment complete: ${mafft_output} \e[0m"
    fi

    #################################################
    # PHASE 5: ANNOTATIONS (PHAROKKA, PHANOTATE, EGG)#
    #################################################

for sample_id in 22 23 24; do
    sample_name="WS2762512A${sample_id}"
    mkdir -p "${ANNOTATION_NITHESIS}/${sample_name}.pharokka"
    mkdir -p "${ANNOTATION_NITHESIS}/${sample_name}.eggnog"

    PHAROKKA_OUT="${ANNOTATION_NITHESIS}/${sample_name}.pharokka"
    EGGNOG_OUT="${ANNOTATION_NITHESIS}/${sample_name}.eggnog"

    # Identify primary contig for annotation
    assembly_nithesis=""
    for cand in "${PHAGETERM_REF}/${sample_name}"_NODE_*.metaSPAdes.fasta "${PHAGETERM_REF}/A${sample_id}"_NODE_*.metaSPAdes.fasta "${SAMPLE_PATH}/phages/assembly_phages/${sample_name}.contigs.filtered.fasta"; do
        if [ -s "${cand}" ]; then
            assembly_nithesis="${cand}"
            break
        fi
    done

    if [ -z "${assembly_nithesis}" ]; then
        echo -e "\e[31m   [${sample_name}] No assembly FASTA found -> skipping annotations \e[0m"
        log_failure "Phase5_Annotations" "${sample_name}" "INPUT_MISSING" "Assembly FASTA not found"
        continue
    fi

    echo -e "\e[31m ====================== \e[0m"
    echo -e "\e[31m PHAROKKA: ${sample_name} \e[0m"
    echo -e "\e[31m ====================== \e[0m"

    expected_pharokka_faa="${PHAROKKA_OUT}/prodigal.faa"
    if [ -s "${expected_pharokka_faa}" ]; then
        echo -e "\e[32m   [${sample_name}] Pharokka output already exists -> skipping \e[0m"
        log_failure "Phase5_Pharokka" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_pharokka_faa} already exists"
    else
        if ! conda run -n BPannotation pharokka.py \
            -t "${threads}" -f \
            -d "${pharokka_db}" \
            -i "${assembly_nithesis}" \
            -g prodigal \
            -o "${PHAROKKA_OUT}" \
            -p "${sample_name}.pharokka"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Pharokka failed \e[0m"
            log_failure "Phase5_Pharokka" "${sample_name}" "EXECUTION_FAILED" "pharokka.py non-zero exit status"
        else
            echo -e "\e[32m Pharokka complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m PHANOTATE: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    expected_phanotate="${ANNOTATION_NITHESIS}/${sample_name}.phanotate"
    if [ -s "${expected_phanotate}" ]; then
        echo -e "\e[32m   [${sample_name}] Phanotate output already exists -> skipping \e[0m"
        log_failure "Phase5_Phanotate" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_phanotate} already exists"
    else
        if ! conda run -n BPannotation phanotate.py \
            --format tabular \
            -o "${expected_phanotate}" \
            "${assembly_nithesis}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Phanotate failed \e[0m"
            log_failure "Phase5_Phanotate" "${sample_name}" "EXECUTION_FAILED" "phanotate.py non-zero exit status"
        else
            echo -e "\e[32m Phanotate complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m =========================== \e[0m"
    echo -e "\e[31m EGGNOG-MAPPER: ${sample_name} \e[0m"
    echo -e "\e[31m =========================== \e[0m"

    expected_eggnog="${EGGNOG_OUT}/${sample_name}.eggnog.emapper.annotations"
    if [ -s "${expected_eggnog}" ]; then
        echo -e "\e[32m   [${sample_name}] EggNOG-mapper output already exists -> skipping \e[0m"
        log_failure "Phase5_EggNOG" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_eggnog} already exists"
    elif [ ! -s "${expected_pharokka_faa}" ]; then
        echo -e "\e[31m   [${sample_name}] Pharokka protein FAA missing -> skipping EggNOG \e[0m"
        log_failure "Phase5_EggNOG" "${sample_name}" "INPUT_MISSING" "Pharokka prodigal.faa missing"
    else
        conda run -n BPannotation emapper.py \
            --cpu "${threads}" \
            --dbmem \
            --data_dir "${eggnog_db}" \
            --temp_dir "${EGGNOG_OUT}" \
            --output_dir "${EGGNOG_OUT}" \
            --no_annot \
            --override \
            -i "${expected_pharokka_faa}" \
            -o "${sample_name}.eggnog" 2>/dev/null || true

        if [ -f "${EGGNOG_OUT}/${sample_name}.eggnog.emapper.seed_orthologs" ]; then
            conda run -n BPannotation emapper.py \
                --cpu "${threads}" \
                --data_dir "${eggnog_db}" \
                --temp_dir "${EGGNOG_OUT}" \
                --output_dir "${EGGNOG_OUT}" \
                --annotate_hits_table "${EGGNOG_OUT}/${sample_name}.eggnog.emapper.seed_orthologs" \
                -m no_search \
                --excel \
                --tax_scope Viruses \
                --override \
                -o "${sample_name}.eggnog" 2>/dev/null || true
            echo -e "\e[32m EggNOG-mapper completed for ${sample_name} \e[0m"
        fi
    fi
done

    #################################################
    # PIPELINE EXECUTION SUMMARY                    #
    #################################################

    echo ""
    echo -e "\e[32m ================================================================ \e[0m"
    echo -e "\e[32m NITHESIS PHAGE ANNOTATION COMPLETE â€” $(date)                     \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
    echo ""
    echo -e "\e[32m -- OUTPUT DIRECTORIES ------------------------------------------ \e[0m"
    echo -e "\e[32m  Genome completeness : ${STRUCTURE_NITHESIS}/*.checkv \e[0m"
    echo -e "\e[32m  Termini & structure : ${STRUCTURE_NITHESIS}/*.term \e[0m"
    echo -e "\e[32m  Annotations         : ${ANNOTATION_NITHESIS}/ \e[0m"
    echo -e "\e[32m  MSA alignment       : ${ANNOTATION_NITHESIS}/mafft_output_nithesis.aln \e[0m"
    echo ""
    echo -e "\e[32m -- LOG FILES --------------------------------------------------- \e[0m"
    echo -e "\e[32m  Full execution log  : ${LOG} \e[0m"

    fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING' "${FAIL_LOG}" 2>/dev/null || echo 0)
    skip_count=$(grep -c '\[SKIPPED' "${FAIL_LOG}" 2>/dev/null || echo 0)

    if [ "${fail_count}" -gt 0 ]; then
        echo -e "\e[31m  Failure/Issues log  : ${FAIL_LOG} (${fail_count} failures detected!) \e[0m"
        echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which samples failed and why. \e[0m"
    else
        echo -e "\e[32m  Failure/Issues log  : ${FAIL_LOG} (0 errors recorded) \e[0m"
    fi
    echo -e "\e[32m  Skipped checkpoints : ${skip_count} records \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
#!/bin/bash
# script for mapping, assembly, and annotations of bacteriophage ILLUMINA samples
# Building No. 1
set -euo pipefail

################
# GLOBAL SETUP #
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/phages"
threads=16

# Paths to databases and host reference
E72_reference="${SAMPLE_PATH}/Escherichia_coli/assembly_Ecoli/WS2762512A08.contigs.filtered.fasta"
checkvdb="${REF_PATH}/checkv-db-v1.5"
eggnog_db="${REF_PATH}/eggnog_db"
pharokka_db="${REF_PATH}/pharokka_db"
truseq2="${REF_PATH}/TruSeq2-PE.fa"
truseq3="${REF_PATH}/TruSeq3-PE.fa"

# Create working directories
PREPROCESSING_PHAGES="${WORK_PATH}/preprocessing_phages"
MAPPING_PHAGES="${WORK_PATH}/mapping_phages"
ASSEMBLY_PHAGES="${WORK_PATH}/assembly_phages"
QUAST_RAW_PHAGES_INPUTS="${ASSEMBLY_PHAGES}/quast_raw_phages_inputs"
STRUCTURE_PHAGES="${WORK_PATH}/structure_phages"
PHAGETERM="${STRUCTURE_PHAGES}/phageterm"
PHAGETERM_REF="${PHAGETERM}/phageterm_references"
ANNOTATION_PHAGES="${WORK_PATH}/annotation_phages"

mkdir -p "${WORK_PATH}"
mkdir -p "${PREPROCESSING_PHAGES}/fastqc_phages_raw"
mkdir -p "${PREPROCESSING_PHAGES}/fastqc_phages_pp"
mkdir -p "${MAPPING_PHAGES}/fastqc_phages_f12"
mkdir -p "${ASSEMBLY_PHAGES}"
mkdir -p "${QUAST_RAW_PHAGES_INPUTS}"
mkdir -p "${STRUCTURE_PHAGES}"
mkdir -p "${PHAGETERM}"
mkdir -p "${PHAGETERM_REF}"
mkdir -p "${ANNOTATION_PHAGES}"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/phages_assembly_mapping_annotation.log"
FAIL_LOG="${WORK_PATH}/phages_assembly_mapping_annotation.failed_skipped.log"
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
echo " Phages assembly, mapping, and annotations - Started: $(date)"
echo " Host reference         : ${E72_reference}"
echo " Full log               : ${LOG}"
echo " Failure & skip log     : ${FAIL_LOG}"
echo "============================================================"

    #################################
    # PHASE 1: PREPROCESSING & QC   #
    #################################

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    read1="${SAMPLE_PATH}/${sample_name}_R1.fastq.gz"
    read2="${SAMPLE_PATH}/${sample_name}_R2.fastq.gz"
    trim_output_1="${PREPROCESSING_PHAGES}/${sample_name}.R1"
    trim_output_2="${PREPROCESSING_PHAGES}/${sample_name}.R2"
    read1t="${trim_output_1}.paired.fastq.gz"
    read2t="${trim_output_2}.paired.fastq.gz"

    if [ ! -s "${read1}" ] || [ ! -s "${read2}" ]; then
        echo -e "\e[31m   [${sample_name}] ERROR: Raw FASTQ missing (${read1} or ${read2}) -> skipping preprocessing \e[0m"
        log_failure "Phase1_Preprocessing" "${sample_name}" "INPUT_MISSING" "Raw fastq missing for ${sample_name}"
        continue
    fi

    echo -e "\e[31m ========================== \e[0m"
    echo -e "\e[31m FASTQC: ${sample_name} (RAW) \e[0m"
    echo -e "\e[31m ========================== \e[0m"

    fastqc_raw_out="${PREPROCESSING_PHAGES}/fastqc_phages_raw/${sample_name}_R1_fastqc.html"
    if [ -s "${fastqc_raw_out}" ]; then
        echo -e "\e[32m   [${sample_name}] FastQC raw already exists -> skipping \e[0m"
        log_failure "Phase1_FastQC_raw" "${sample_name}" "SKIPPED_EXISTS" "Output ${fastqc_raw_out} already exists"
    else
        conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${PREPROCESSING_PHAGES}/fastqc_phages_raw" \
            "${read1}" "${read2}" 2>/dev/null || true
    fi

    echo -e "\e[31m ========================= \e[0m"
    echo -e "\e[31m TRIMMOMATIC: ${sample_name} \e[0m"
    echo -e "\e[31m ========================= \e[0m"

    if [ -s "${read1t}" ] && [ -s "${read2t}" ]; then
        echo -e "\e[32m   [${sample_name}] Trimmed FASTQ already exists -> skipping Trimmomatic \e[0m"
        log_failure "Phase1_Trimmomatic" "${sample_name}" "SKIPPED_EXISTS" "Output ${read1t} and ${read2t} already exist"
    else
        if ! conda run -n preprocessing trimmomatic PE \
            -threads "${threads}" \
            -phred64 \
            "${read1}" "${read2}" \
            "${trim_output_1}.paired.fastq.gz" "${trim_output_1}.unpaired.fastq.gz" \
            "${trim_output_2}.paired.fastq.gz" "${trim_output_2}.unpaired.fastq.gz" \
            ILLUMINACLIP:"${truseq2}":2:30:10:2:True \
            ILLUMINACLIP:"${truseq3}":2:30:10:2:True \
            LEADING:3 TRAILING:3 MINLEN:50; then
            echo -e "\e[31m   [${sample_name}] ERROR: Trimmomatic execution failed \e[0m"
            log_failure "Phase1_Trimmomatic" "${sample_name}" "EXECUTION_FAILED" "trimmomatic non-zero exit status"
        else
            echo -e "\e[32m Trimmomatic complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ========================= \e[0m"
    echo -e "\e[31m FASTQC: ${sample_name} (PP) \e[0m"
    echo -e "\e[31m ========================= \e[0m"

    fastqc_pp_out="${PREPROCESSING_PHAGES}/fastqc_phages_pp/${sample_name}.R1.paired_fastqc.html"
    if [ -s "${fastqc_pp_out}" ]; then
        echo -e "\e[32m   [${sample_name}] FastQC pp already exists -> skipping \e[0m"
        log_failure "Phase1_FastQC_pp" "${sample_name}" "SKIPPED_EXISTS" "Output ${fastqc_pp_out} already exists"
    elif [ -s "${read1t}" ] && [ -s "${read2t}" ]; then
        conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${PREPROCESSING_PHAGES}/fastqc_phages_pp" \
            "${read1t}" "${read2t}" 2>/dev/null || true
    fi
done

    echo -e "\e[31m ============================ \e[0m"
    echo -e "\e[31m MULTIQC: PHAGE SAMPLES (RAW) \e[0m"
    echo -e "\e[31m ============================ \e[0m"

    conda run -n preprocessing multiqc \
        --force \
        "${PREPROCESSING_PHAGES}/fastqc_phages_raw" \
        --filename "${PREPROCESSING_PHAGES}/multiqc_phages_raw" 2>/dev/null || true

    echo -e "\e[31m =========================== \e[0m"
    echo -e "\e[31m MULTIQC: PHAGE SAMPLES (PP) \e[0m"
    echo -e "\e[31m =========================== \e[0m"

    conda run -n preprocessing multiqc \
        --force \
        "${PREPROCESSING_PHAGES}/fastqc_phages_pp" \
        --filename "${PREPROCESSING_PHAGES}/multiqc_phages_pp" 2>/dev/null || true

    #################################################
    # PHASE 2: HOST READ MAPPING & FILTERING (-f 12)#
    #################################################

    if [ -s "${E72_reference}" ]; then
        if [ ! -f "${E72_reference}.bwt" ]; then
            echo -e "\e[32m Indexing host reference: ${E72_reference} \e[0m"
            conda run -n mapping bwa index "${E72_reference}"
        fi
    else
        echo -e "\e[33m WARNING: Host reference not found: ${E72_reference} \e[0m"
        log_failure "Phase2_HostMapping" "HOST" "INPUT_MISSING" "Host reference missing: ${E72_reference}"
    fi

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    MAPPING_PHAGE_SAMPLE="${MAPPING_PHAGES}/${sample_name}"
    mkdir -p "${MAPPING_PHAGE_SAMPLE}"

    read1t="${PREPROCESSING_PHAGES}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_PHAGES}/${sample_name}.R2.paired.fastq.gz"
    mapped_bam="${MAPPING_PHAGE_SAMPLE}/${sample_name}.bam"
    sorted_bam="${MAPPING_PHAGE_SAMPLE}/${sample_name}.sorted_coord.bam"
    unmapped_bam="${MAPPING_PHAGE_SAMPLE}/${sample_name}.f12.bam"
    unmapped_name_bam="${MAPPING_PHAGE_SAMPLE}/${sample_name}.sorted_name.f12.bam"
    read1f="${MAPPING_PHAGE_SAMPLE}/${sample_name}.R1.f12.fastq.gz"
    read2f="${MAPPING_PHAGE_SAMPLE}/${sample_name}.R2.f12.fastq.gz"

    if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        echo -e "\e[31m   [${sample_name}] Trimmed reads missing -> skipping host subtraction \e[0m"
        log_failure "Phase2_HostMapping" "${sample_name}" "INPUT_MISSING" "Trimmed reads missing"
        continue
    fi

    if [ -s "${read1f}" ] && [ -s "${read2f}" ]; then
        echo -e "\e[32m   [${sample_name}] Unmapped phage FASTQs already generated -> skipping BWA host mapping \e[0m"
        log_failure "Phase2_HostMapping" "${sample_name}" "SKIPPED_EXISTS" "Output ${read1f} and ${read2f} already exist"
    elif [ ! -s "${E72_reference}" ]; then
        echo -e "\e[33m   [${sample_name}] Host reference missing, using trimmed reads directly \e[0m"
        cp -f "${read1t}" "${read1f}"
        cp -f "${read2t}" "${read2f}"
    else
        echo -e "\e[31m ==================== \e[0m"
        echo -e "\e[31m BWA MEM: ${sample_name} \e[0m"
        echo -e "\e[31m ==================== \e[0m"

        conda run -n mapping bwa mem \
            -t "${threads}" \
            "${E72_reference}" \
            "${read1t}" "${read2t}" \
            > "${MAPPING_PHAGE_SAMPLE}/${sample_name}.sam"

        conda run -n mapping samtools view \
            -@ "${threads}" \
            -bS \
            -o "${mapped_bam}" \
            "${MAPPING_PHAGE_SAMPLE}/${sample_name}.sam"

        rm -f "${MAPPING_PHAGE_SAMPLE}/${sample_name}.sam"

        conda run -n mapping samtools sort \
            -@ "${threads}" \
            -o "${sorted_bam}" \
            "${mapped_bam}"

        conda run -n mapping samtools index \
            -@ "${threads}" \
            "${sorted_bam}"

        conda run -n mapping samtools coverage \
            "${sorted_bam}" \
            > "${MAPPING_PHAGE_SAMPLE}/${sample_name}.sorted_coord.coverage.txt"

        # Filter unmapped reads (flag 12 = both reads of pair unmapped)
        conda run -n mapping samtools view \
            -@ "${threads}" \
            -f 12 \
            -b \
            "${sorted_bam}" \
            -o "${unmapped_bam}"

        conda run -n mapping samtools sort \
            -@ "${threads}" \
            -n \
            -o "${unmapped_name_bam}" \
            "${unmapped_bam}"

        conda run -n mapping samtools bam2fq \
            -@ "${threads}" -n \
            -1 "${read1f}" \
            -2 "${read2f}" \
            "${unmapped_name_bam}"

        echo -e "\e[32m Unmapped phage reads extracted: ${read1f}, ${read2f} \e[0m"
    fi

    # QC on filtered phage reads
    fastqc_f12_out="${MAPPING_PHAGES}/fastqc_phages_f12/${sample_name}.R1.f12_fastqc.html"
    if [ -s "${fastqc_f12_out}" ]; then
        echo -e "\e[32m   [${sample_name}] FastQC f12 already exists -> skipping \e[0m"
        log_failure "Phase2_FastQC_f12" "${sample_name}" "SKIPPED_EXISTS" "Output ${fastqc_f12_out} already exists"
    elif [ -s "${read1f}" ] && [ -s "${read2f}" ]; then
        conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${MAPPING_PHAGES}/fastqc_phages_f12" \
            "${read1f}" "${read2f}" 2>/dev/null || true
    fi
done

    echo -e "\e[31m ============================== \e[0m"
    echo -e "\e[31m MULTIQC: PHAGE SAMPLES (-f 12) \e[0m"
    echo -e "\e[31m ============================== \e[0m"

    conda run -n preprocessing multiqc \
        --force \
        "${MAPPING_PHAGES}/fastqc_phages_f12" \
        --filename "${MAPPING_PHAGES}/multiqc_phages_f12" 2>/dev/null || true

    #################################################
    # PHASE 3: DE NOVO ASSEMBLY & QC (SPADES META)  #
    #################################################

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    MAPPING_PHAGE_SAMPLE="${MAPPING_PHAGES}/${sample_name}"
    read1f="${MAPPING_PHAGE_SAMPLE}/${sample_name}.R1.f12.fastq.gz"
    read2f="${MAPPING_PHAGE_SAMPLE}/${sample_name}.R2.f12.fastq.gz"
    assembly="${ASSEMBLY_PHAGES}/${sample_name}.contig.fasta/contigs.fasta"
    filtered_assembly="${ASSEMBLY_PHAGES}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${read1f}" ] || [ ! -s "${read2f}" ]; then
        echo -e "\e[31m   [${sample_name}] Phage reads missing -> skipping de novo assembly \e[0m"
        log_failure "Phase3_Assembly" "${sample_name}" "INPUT_MISSING" "Filtered phage reads missing"
        continue
    fi

    # SeqKit coverage check on input reads
    conda run -n assembly seqkit stats \
        --threads "${threads}" \
        --all --tabular \
        -o "${MAPPING_PHAGE_SAMPLE}/${sample_name}.f12.fastq.stats.txt" \
        "${read1f}" "${read2f}" 2>/dev/null || true

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m SPADES: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${assembly}" ]; then
        echo -e "\e[32m   [${sample_name}] metaSPAdes assembly already exists -> skipping \e[0m"
        log_failure "Phase3_metaSPAdes" "${sample_name}" "SKIPPED_EXISTS" "Output ${assembly} already exists"
    else
        if ! conda run -n assembly spades.py \
            --meta \
            -1 "${read1f}" \
            -2 "${read2f}" \
            -t "${threads}" \
            -o "${ASSEMBLY_PHAGES}/${sample_name}.contig.fasta"; then
            echo -e "\e[31m   [${sample_name}] ERROR: metaSPAdes failed \e[0m"
            log_failure "Phase3_metaSPAdes" "${sample_name}" "EXECUTION_FAILED" "spades.py non-zero exit status"
            continue
        else
            echo -e "\e[32m metaSPAdes complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ========================= \e[0m"
    echo -e "\e[31m QUAST: ${sample_name} (RAW) \e[0m"
    echo -e "\e[31m ========================= \e[0m"

    if [ -s "${assembly}" ]; then
        conda run -n assembly quast.py \
            "${assembly}" \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_PHAGES}/${sample_name}.raw.quast" 2>/dev/null || true
    fi

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m SEQKIT: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${filtered_assembly}" ]; then
        echo -e "\e[32m   [${sample_name}] Filtered contigs already exist -> skipping SeqKit \e[0m"
        log_failure "Phase3_SeqKit" "${sample_name}" "SKIPPED_EXISTS" "Output ${filtered_assembly} already exists"
    elif [ -s "${assembly}" ]; then
        conda run -n assembly seqkit seq \
            --min-len 500 \
            "${assembly}" \
            > "${filtered_assembly}"
    fi

    echo -e "\e[31m ======================== \e[0m"
    echo -e "\e[31m QUAST: ${sample_name} (PP) \e[0m"
    echo -e "\e[31m ======================== \e[0m"

    if [ -s "${filtered_assembly}" ]; then
        conda run -n assembly quast.py \
            "${filtered_assembly}" \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_PHAGES}/${sample_name}.filtered.quast" 2>/dev/null || true
    fi
done

    echo -e "\e[31m ============================== \e[0m"
    echo -e "\e[31m QUAST: ALL PHAGE SAMPLES (RAW) \e[0m"
    echo -e "\e[31m ============================== \e[0m"

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    if [ -s "${ASSEMBLY_PHAGES}/${sample_name}.contig.fasta/contigs.fasta" ]; then
        cp -f "${ASSEMBLY_PHAGES}/${sample_name}.contig.fasta/contigs.fasta" \
            "${QUAST_RAW_PHAGES_INPUTS}/${sample_name}.contigs.raw.fasta"
    fi
done

    if compgen -G "${QUAST_RAW_PHAGES_INPUTS}/*.fasta" > /dev/null; then
        conda run -n assembly quast.py \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_PHAGES}/all_phages.raw.quast" \
            "${QUAST_RAW_PHAGES_INPUTS}/"*.contigs.raw.fasta 2>/dev/null || true
    fi

    echo -e "\e[31m ============================= \e[0m"
    echo -e "\e[31m QUAST: ALL PHAGE SAMPLES (PP) \e[0m"
    echo -e "\e[31m ============================= \e[0m"

    if compgen -G "${ASSEMBLY_PHAGES}/*.contigs.filtered.fasta" > /dev/null; then
        conda run -n assembly quast.py \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_PHAGES}/all_phages.filtered.quast" \
            "${ASSEMBLY_PHAGES}/"*.contigs.filtered.fasta 2>/dev/null || true
    fi

    #################################################
    # PHASE 4: CHECK GENOME COMPLETENESS (CHECKV)   #
    #################################################

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    filtered_assembly="${ASSEMBLY_PHAGES}/${sample_name}.contigs.filtered.fasta"
    checkv_dir="${STRUCTURE_PHAGES}/${sample_name}.checkv"
    checkv_summary="${checkv_dir}/quality_summary.tsv"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m CHECKV: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${checkv_summary}" ]; then
        echo -e "\e[32m   [${sample_name}] CheckV output already exists -> skipping \e[0m"
        log_failure "Phase4_CheckV" "${sample_name}" "SKIPPED_EXISTS" "Output ${checkv_summary} already exists"
    elif [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Filtered assembly missing -> skipping CheckV \e[0m"
        log_failure "Phase4_CheckV" "${sample_name}" "INPUT_MISSING" "Filtered assembly missing"
    else
        rm -rf "${checkv_dir}"
        if ! conda run -n BPstructure checkv end_to_end \
            -t "${threads}" \
            -d "${checkvdb}" \
            "${filtered_assembly}" \
            "${checkv_dir}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: CheckV failed \e[0m"
            log_failure "Phase4_CheckV" "${sample_name}" "EXECUTION_FAILED" "checkv non-zero exit status"
        else
            echo -e "\e[32m CheckV complete for ${sample_name} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 5: PHAGE TERMINI PREDICTION (PHAGETERM) #
    #################################################

    echo -e "\e[32m Extracting high-completeness contigs (>=90%) from CheckV results... \e[0m"

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    checkv_tsv="${STRUCTURE_PHAGES}/${sample_name}.checkv/quality_summary.tsv"
    filtered_assembly="${ASSEMBLY_PHAGES}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${checkv_tsv}" ] || [ ! -s "${filtered_assembly}" ]; then
        continue
    fi

    tail -n +2 "${checkv_tsv}" \
    | awk -F'\t' '($9 != "NA" && $9+0 >= 90) || ($10 != "NA" && $10+0 >= 90) { print $1 }' \
    | while read -r contig_id; do
        node_num=$(echo "${contig_id}" | grep -oP '(?<=NODE_)\d+' || true)
        if [ -z "${node_num}" ]; then
            continue
        fi

        out_fasta="${PHAGETERM_REF}/${sample_name}_NODE_${node_num}.metaSPAdes.fasta"
        new_header="${sample_name}_NODE_${node_num}"

        if [ -s "${out_fasta}" ]; then
            echo -e "\e[32m   Node ${node_num} already extracted for ${sample_name} -> skipping \e[0m"
            log_failure "Phase5_NodeExtract" "${sample_name}" "SKIPPED_EXISTS" "File ${out_fasta} already exists"
        else
            tmp_fa="${PHAGETERM_REF}/tmp.${sample_name}.${node_num}.fa"
            conda run -n assembly seqkit grep -p "${contig_id}" "${filtered_assembly}" > "${tmp_fa}"
            conda run -n assembly seqkit replace -p '^.*$' -r "${new_header}" "${tmp_fa}" > "${out_fasta}"
            rm -f "${tmp_fa}"
            echo -e "\e[32m   Saved: ${out_fasta} \e[0m"
        fi
    done
done

for phageterm_ref in "${PHAGETERM_REF}"/*.metaSPAdes.fasta; do
    [ -f "${phageterm_ref}" ] || continue

    basename_no_ext=$(basename "${phageterm_ref}" .metaSPAdes.fasta)
    sample_name=$(echo "${basename_no_ext}" | grep -oP '^WS\d+A\d+' || true)
    node_num=$(echo "${basename_no_ext}" | grep -oP '(?<=NODE_)\d+' || true)

    if [ -z "${sample_name}" ] || [ -z "${node_num}" ]; then
        continue
    fi

    read1f="${MAPPING_PHAGES}/${sample_name}/${sample_name}.R1.f12.fastq.gz"
    read2f="${MAPPING_PHAGES}/${sample_name}/${sample_name}.R2.f12.fastq.gz"

    if [ ! -s "${read1f}" ] || [ ! -s "${read2f}" ]; then
        echo -e "\e[31m   [${sample_name}_NODE_${node_num}] Reads missing -> skipping PhageTerm \e[0m"
        log_failure "Phase5_PhageTerm" "${sample_name}_NODE_${node_num}" "INPUT_MISSING" "Filtered reads missing"
        continue
    fi

    TERM_PHAGE="${STRUCTURE_PHAGES}/${sample_name}_NODE_${node_num}.term"
    mkdir -p "${TERM_PHAGE}"

    echo -e "\e[31m ================================ \e[0m"
    echo -e "\e[31m PHAGETERM: ${sample_name} NODE_${node_num} \e[0m"
    echo -e "\e[31m ================================ \e[0m"

    expected_term="${TERM_PHAGE}/${sample_name}_NODE_${node_num}_report.pdf"
    if [ -s "${expected_term}" ]; then
        echo -e "\e[32m   [${sample_name}_NODE_${node_num}] PhageTerm output already exists -> skipping \e[0m"
        log_failure "Phase5_PhageTerm" "${sample_name}_NODE_${node_num}" "SKIPPED_EXISTS" "Output ${expected_term} already exists"
    else
        if ! conda run -n phageterm_env phageterm \
            -f "${read1f}" \
            -p "${read2f}" \
            -r "${phageterm_ref}" \
            -c "${threads}" \
            -o "${TERM_PHAGE}" \
            --report_title "${sample_name}_NODE_${node_num}"; then
            echo -e "\e[31m   [${sample_name}_NODE_${node_num}] ERROR: PhageTerm failed \e[0m"
            log_failure "Phase5_PhageTerm" "${sample_name}_NODE_${node_num}" "EXECUTION_FAILED" "phageterm non-zero exit status"
        else
            echo -e "\e[32m PhageTerm complete: ${sample_name} NODE_${node_num} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 6: PHAGE ANNOTATIONS (PHAROKKA, EGGNOG) #
    #################################################

for sample_id in {21..24}; do
    sample_name="WS2762512A${sample_id}"
    filtered_assembly="${WORK_PATH}/assembly_phages/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Filtered assembly missing -> skipping annotations \e[0m"
        log_failure "Phase6_Annotations" "${sample_name}" "INPUT_MISSING" "Filtered assembly missing"
        continue
    fi

    mkdir -p "${ANNOTATION_PHAGES}/${sample_name}.pharokka"
    mkdir -p "${ANNOTATION_PHAGES}/${sample_name}.eggnog"

    PHAROKKA_OUT="${ANNOTATION_PHAGES}/${sample_name}.pharokka"
    EGGNOG_PHAGE="${ANNOTATION_PHAGES}/${sample_name}.eggnog"
    prodigal_input="${PHAROKKA_OUT}/prodigal.faa"

    echo -e "\e[31m ====================== \e[0m"
    echo -e "\e[31m PHAROKKA: ${sample_name} \e[0m"
    echo -e "\e[31m ====================== \e[0m"

    if [ -s "${prodigal_input}" ]; then
        echo -e "\e[32m   [${sample_name}] Pharokka annotation already exists -> skipping \e[0m"
        log_failure "Phase6_Pharokka" "${sample_name}" "SKIPPED_EXISTS" "Output ${prodigal_input} already exists"
    else
        if ! conda run -n BPannotation pharokka.py \
            -t "${threads}" -f \
            -d "${pharokka_db}" \
            -i "${filtered_assembly}" \
            -g prodigal \
            -o "${PHAROKKA_OUT}" \
            -p "${sample_name}.pharokka"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Pharokka failed \e[0m"
            log_failure "Phase6_Pharokka" "${sample_name}" "EXECUTION_FAILED" "pharokka.py non-zero exit status"
        else
            echo -e "\e[32m Pharokka complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m PHANOTATE: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    expected_phanotate="${ANNOTATION_PHAGES}/${sample_name}.phanotate"
    if [ -s "${expected_phanotate}" ]; then
        echo -e "\e[32m   [${sample_name}] Phanotate already exists -> skipping \e[0m"
        log_failure "Phase6_Phanotate" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_phanotate} already exists"
    else
        if ! conda run -n BPannotation phanotate.py \
            --format tabular \
            -o "${expected_phanotate}" \
            "${filtered_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Phanotate failed \e[0m"
            log_failure "Phase6_Phanotate" "${sample_name}" "EXECUTION_FAILED" "phanotate.py non-zero exit status"
        else
            echo -e "\e[32m Phanotate complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m =========================== \e[0m"
    echo -e "\e[31m EGGNOG-MAPPER: ${sample_name} \e[0m"
    echo -e "\e[31m =========================== \e[0m"

    expected_eggnog="${EGGNOG_PHAGE}/${sample_name}.eggnog.emapper.annotations"
    if [ -s "${expected_eggnog}" ]; then
        echo -e "\e[32m   [${sample_name}] EggNOG-mapper already exists -> skipping \e[0m"
        log_failure "Phase6_EggNOG" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_eggnog} already exists"
    elif [ ! -s "${prodigal_input}" ]; then
        echo -e "\e[31m   [${sample_name}] Pharokka prodigal.faa missing -> skipping EggNOG \e[0m"
        log_failure "Phase6_EggNOG" "${sample_name}" "INPUT_MISSING" "Pharokka prodigal.faa missing"
    else
        conda run -n BPannotation emapper.py \
            --cpu "${threads}" \
            --dbmem \
            --data_dir "${eggnog_db}" \
            --temp_dir "${EGGNOG_PHAGE}" \
            --output_dir "${EGGNOG_PHAGE}" \
            --no_annot \
            --override \
            -i "${prodigal_input}" \
            -o "${sample_name}.eggnog" 2>/dev/null || true

        if [ -f "${EGGNOG_PHAGE}/${sample_name}.eggnog.emapper.seed_orthologs" ]; then
            conda run -n BPannotation emapper.py \
                --cpu "${threads}" \
                --data_dir "${eggnog_db}" \
                --temp_dir "${EGGNOG_PHAGE}" \
                --output_dir "${EGGNOG_PHAGE}" \
                --annotate_hits_table "${EGGNOG_PHAGE}/${sample_name}.eggnog.emapper.seed_orthologs" \
                -m no_search \
                --excel \
                --tax_scope Viruses \
                --override \
                -o "${sample_name}.eggnog" 2>/dev/null || true
            echo -e "\e[32m EggNOG-mapper completed for sample ${sample_name} \e[0m"
        fi
    fi
done

    #################################################
    # PIPELINE EXECUTION SUMMARY                    #
    #################################################

    echo ""
    echo -e "\e[32m ================================================================ \e[0m"
    echo -e "\e[32m PHAGES ASSEMBLY, MAPPING & ANNOTATION COMPLETE â€” $(date)         \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
    echo ""
    echo -e "\e[32m -- OUTPUT DIRECTORIES ------------------------------------------ \e[0m"
    echo -e "\e[32m  Preprocessing       : ${PREPROCESSING_PHAGES} \e[0m"
    echo -e "\e[32m  Mapping             : ${MAPPING_PHAGES} \e[0m"
    echo -e "\e[32m  Phage assemblies    : ${ASSEMBLY_PHAGES} \e[0m"
    echo -e "\e[32m  Annotation          : ${ANNOTATION_PHAGES} \e[0m"
    echo -e "\e[32m  Genome completeness : ${STRUCTURE_PHAGES}/*.checkv \e[0m"
    echo -e "\e[32m  Termini & structure : ${STRUCTURE_PHAGES}/*.term \e[0m"
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
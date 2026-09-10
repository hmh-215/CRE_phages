#!/bin/bash
# scripts for de novo assembly, annotation, typing, AMR and virulence profiling of K. pneumoniae ILLUMINA samples
# Building No. 1 (Project 2)
set -euo pipefail

################
# GLOBAL SETUP #
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages_2"
WORK_PATH="${SAMPLE_PATH}/Klebsiella_pneumoniae"

threads=16

# Paths to databases
truseq2="${REF_PATH}/TruSeq2-PE.fa"
truseq3="${REF_PATH}/TruSeq3-PE.fa"
bakta_db="${REF_PATH}/bakta_db/db"
amrfinder_db="${REF_PATH}/bakta_db/db/amrfinderplus-db/latest"
eggnog_db="${REF_PATH}/eggnog_db"

# Create working directories
PREPROCESSING_KLEB="${WORK_PATH}/preprocessing_Kleb"
ASSEMBLY_KLEB="${WORK_PATH}/assembly_Kleb"
ANNOTATION_KLEB="${WORK_PATH}/annotation_Kleb"
ANNOTATION_KLEB_CHECKM="${ANNOTATION_KLEB}/checkm"
ANNOTATION_KLEB_MLST="${ANNOTATION_KLEB}/mlst"
QUAST_RAW_KLEB_INPUTS="${ASSEMBLY_KLEB}/quast_raw_Kleb_inputs"
CHECKM_INPUTS="${ANNOTATION_KLEB}/checkm/inputs_Kleb"
ANNOTATION_KLEB_BAKTA_BASE="${ANNOTATION_KLEB}/bakta"
ANNOTATION_KLEB_AMRFINDER_BASE="${ANNOTATION_KLEB}/amrfinder"
ANNOTATION_KLEB_EGGNOG_BASE="${ANNOTATION_KLEB}/eggnog"
AMR_KLEB="${WORK_PATH}/amr_Kleb"
ABRICATE_KLEB="${AMR_KLEB}/abricate"
PAI_KLEB="${WORK_PATH}/pai_Kleb"
KLEBORATE_KLEB="${PAI_KLEB}/kleborate"

mkdir -p "${WORK_PATH}"
mkdir -p "${PREPROCESSING_KLEB}/fastqc_Kleb_raw"
mkdir -p "${PREPROCESSING_KLEB}/fastqc_Kleb_pp"
mkdir -p "${ASSEMBLY_KLEB}"
mkdir -p "${ANNOTATION_KLEB}"
mkdir -p "${QUAST_RAW_KLEB_INPUTS}"
mkdir -p "${ANNOTATION_KLEB_CHECKM}"
mkdir -p "${CHECKM_INPUTS}"
mkdir -p "${ANNOTATION_KLEB_MLST}"
mkdir -p "${ANNOTATION_KLEB_BAKTA_BASE}"
mkdir -p "${ANNOTATION_KLEB_AMRFINDER_BASE}"
mkdir -p "${ANNOTATION_KLEB_EGGNOG_BASE}"
mkdir -p "${ABRICATE_KLEB}"
mkdir -p "${KLEBORATE_KLEB}"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/Kleb_assembly_annotation_2.log"
FAIL_LOG="${WORK_PATH}/Kleb_assembly_annotation_2.failed_skipped.log"
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

sample_ids=(K10-LJF30359 K14-LJF30358 K16-LJF30357 K66-LJF30363 K74-LJF30355 K75-LJF30356 EM10311-LJF30362)

echo "============================================================"
echo " K. pneumoniae assemblies and annotations (Project 2) - Started: $(date)"
echo " Total samples cohort   : ${#sample_ids[@]}"
echo " Full log               : ${LOG}"
echo " Failure & skip log     : ${FAIL_LOG}"
echo "============================================================"

    #################################
    # PHASE 1: PREPROCESSING & QC   #
    #################################

for sample_id in "${sample_ids[@]}"; do
    sample_name="${sample_id}"
    read1="${SAMPLE_PATH}/${sample_id}_L3_1.fq.gz"
    read2="${SAMPLE_PATH}/${sample_id}_L3_2.fq.gz"

    trim_output_1="${PREPROCESSING_KLEB}/${sample_id}.R1"
    trim_output_2="${PREPROCESSING_KLEB}/${sample_id}.R2"
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

    fastqc_raw_out="${PREPROCESSING_KLEB}/fastqc_Kleb_raw/${sample_id}_L3_1_fastqc.html"
    if [ -s "${fastqc_raw_out}" ]; then
        echo -e "\e[32m   [${sample_name}] FastQC raw already exists -> skipping \e[0m"
        log_failure "Phase1_FastQC_raw" "${sample_name}" "SKIPPED_EXISTS" "Output ${fastqc_raw_out} already exists"
    else
        conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${PREPROCESSING_KLEB}/fastqc_Kleb_raw" \
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
            -phred33 \
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

    fastqc_pp_out="${PREPROCESSING_KLEB}/fastqc_Kleb_pp/${sample_id}.R1.paired_fastqc.html"
    if [ -s "${fastqc_pp_out}" ]; then
        echo -e "\e[32m   [${sample_name}] FastQC pp already exists -> skipping \e[0m"
        log_failure "Phase1_FastQC_pp" "${sample_name}" "SKIPPED_EXISTS" "Output ${fastqc_pp_out} already exists"
    elif [ -s "${read1t}" ] && [ -s "${read2t}" ]; then
        conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${PREPROCESSING_KLEB}/fastqc_Kleb_pp" \
            "${read1t}" "${read2t}" 2>/dev/null || true
    fi
done

    echo -e "\e[31m ===================================== \e[0m"
    echo -e "\e[31m MULTIQC: ALL K. PNEUMONIAE SAMPLES (RAW) \e[0m"
    echo -e "\e[31m ===================================== \e[0m"

    conda run -n preprocessing multiqc \
        --force \
        "${PREPROCESSING_KLEB}/fastqc_Kleb_raw" \
        --filename "${PREPROCESSING_KLEB}/multiqc_Kleb_raw" 2>/dev/null || true

    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m MULTIQC: ALL K. PNEUMONIAE SAMPLES (PP) \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    conda run -n preprocessing multiqc \
        --force \
        "${PREPROCESSING_KLEB}/fastqc_Kleb_pp" \
        --filename "${PREPROCESSING_KLEB}/multiqc_Kleb_pp" 2>/dev/null || true

    #################################
    # PHASE 2: DE NOVO ASSEMBLY & QC#
    #################################

for sample_id in "${sample_ids[@]}"; do
    sample_name="${sample_id}"
    read1t="${PREPROCESSING_KLEB}/${sample_id}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_KLEB}/${sample_id}.R2.paired.fastq.gz"

    if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        echo -e "\e[31m   [${sample_name}] Trimmed FASTQ missing -> skipping assembly \e[0m"
        log_failure "Phase2_SPAdes" "${sample_name}" "INPUT_MISSING" "Trimmed reads missing for ${sample_name}"
        continue
    fi

    Kleb_assembly="${ASSEMBLY_KLEB}/${sample_name}.contig.fasta/contigs.fasta"
    filtered_Kleb_assembly="${ASSEMBLY_KLEB}/${sample_name}.contigs.filtered.fasta"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m SPADES: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${Kleb_assembly}" ]; then
        echo -e "\e[32m   [${sample_name}] SPAdes assembly already exists -> skipping \e[0m"
        log_failure "Phase2_SPAdes" "${sample_name}" "SKIPPED_EXISTS" "Output ${Kleb_assembly} already exists"
    else
        if ! conda run -n assembly spades.py \
            -1 "${read1t}" \
            -2 "${read2t}" \
            --careful \
            -t "${threads}" \
            -o "${ASSEMBLY_KLEB}/${sample_name}.contig.fasta"; then
            echo -e "\e[31m   [${sample_name}] ERROR: SPAdes failed \e[0m"
            log_failure "Phase2_SPAdes" "${sample_name}" "EXECUTION_FAILED" "spades.py non-zero exit status"
            continue
        else
            echo -e "\e[32m SPAdes complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ========================= \e[0m"
    echo -e "\e[31m QUAST: ${sample_name} (RAW) \e[0m"
    echo -e "\e[31m ========================= \e[0m"

    if [ -s "${Kleb_assembly}" ]; then
        conda run -n assembly quast.py \
            "${Kleb_assembly}" \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_KLEB}/${sample_name}.raw.quast" 2>/dev/null || true
    fi

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m SEQKIT: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${filtered_Kleb_assembly}" ]; then
        echo -e "\e[32m   [${sample_name}] Filtered contigs already exist -> skipping SeqKit \e[0m"
        log_failure "Phase2_SeqKit" "${sample_name}" "SKIPPED_EXISTS" "Output ${filtered_Kleb_assembly} already exists"
    elif [ -s "${Kleb_assembly}" ]; then
        conda run -n assembly seqkit seq \
            --min-len 500 \
            "${Kleb_assembly}" \
            > "${filtered_Kleb_assembly}"
    fi

    echo -e "\e[31m ======================== \e[0m"
    echo -e "\e[31m QUAST: ${sample_name} (PP) \e[0m"
    echo -e "\e[31m ======================== \e[0m"

    if [ -s "${filtered_Kleb_assembly}" ]; then
        conda run -n assembly quast.py \
            "${filtered_Kleb_assembly}" \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_KLEB}/${sample_name}.filtered.quast" 2>/dev/null || true
    fi
done

    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m QUAST: ALL K. PNEUMONIAE SAMPLES (RAW) \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${sample_id}"
    if [ -s "${ASSEMBLY_KLEB}/${sample_name}.contig.fasta/contigs.fasta" ]; then
        cp -f "${ASSEMBLY_KLEB}/${sample_name}.contig.fasta/contigs.fasta" "${QUAST_RAW_KLEB_INPUTS}/${sample_name}.contigs.raw.fasta"
    fi
done

    if compgen -G "${QUAST_RAW_KLEB_INPUTS}/*.fasta" > /dev/null; then
        conda run -n assembly quast.py \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_KLEB}/all_Kleb.raw.quast" \
            "${QUAST_RAW_KLEB_INPUTS}/"*.contigs.raw.fasta 2>/dev/null || true
    fi

    echo -e "\e[31m =================================== \e[0m"
    echo -e "\e[31m QUAST: ALL K. PNEUMONIAE SAMPLES (PP) \e[0m"
    echo -e "\e[31m =================================== \e[0m"

    if compgen -G "${ASSEMBLY_KLEB}/*contigs.filtered.fasta" > /dev/null; then
        conda run -n assembly quast.py \
            --threads "${threads}" \
            --output-dir "${ASSEMBLY_KLEB}/all_Kleb.filtered.quast" \
            "${ASSEMBLY_KLEB}/"*contigs.filtered.fasta 2>/dev/null || true
    fi

    #################################
    # PHASE 3: CHECKM COMPLETENESS  #
    #################################

    for f in "${ASSEMBLY_KLEB}/"*contigs.filtered.fasta; do
        [ -f "${f}" ] || continue
        cp -f "${f}" "${CHECKM_INPUTS}/"
    done

    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m CHECKM: ALL K. PNEUMONIAE SAMPLES (PP) \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    expected_checkm="${ANNOTATION_KLEB_CHECKM}/all_Kleb.quality.checkm.tsv"
    if [ -s "${expected_checkm}" ]; then
        echo -e "\e[32m CheckM summary already exists -> skipping \e[0m"
        log_failure "Phase3_CheckM" "ALL" "SKIPPED_EXISTS" "Output ${expected_checkm} already exists"
    elif compgen -G "${CHECKM_INPUTS}/*.fasta" > /dev/null; then
        conda run -n BPstructure checkm lineage_wf \
            --threads "${threads}" \
            --extension fasta \
            --tab_table \
            "${CHECKM_INPUTS}" \
            "${ANNOTATION_KLEB_CHECKM}" 2>/dev/null || true

        if [ -f "${ANNOTATION_KLEB_CHECKM}/lineage.ms" ]; then
            conda run -n BPstructure checkm qa \
                "${ANNOTATION_KLEB_CHECKM}/lineage.ms" \
                "${ANNOTATION_KLEB_CHECKM}" \
                -o 2 \
                --tab_table \
                -f "${expected_checkm}"
            echo -e "\e[32m CheckM complete: ${expected_checkm} \e[0m"
        else
            echo -e "\e[33m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
        fi
    fi

    #################################
    # PHASE 4: MLST & KLEBORATE     #
    #################################

    echo -e "\e[31m ============================= \e[0m"
    echo -e "\e[31m MLST: ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ============================= \e[0m"

    expected_mlst="${ANNOTATION_KLEB_MLST}/all_Kleb.mlst.csv"
    if [ -s "${expected_mlst}" ]; then
        echo -e "\e[32m MLST output already exists -> skipping \e[0m"
        log_failure "Phase4_MLST" "ALL" "SKIPPED_EXISTS" "Output ${expected_mlst} already exists"
    elif compgen -G "${ASSEMBLY_KLEB}/*contigs.filtered.fasta" > /dev/null; then
        conda run -n BPtyping mlst \
            --threads "${threads}" \
            --scheme klebsiella \
            --csv \
            "${ASSEMBLY_KLEB}/"*contigs.filtered.fasta \
            > "${expected_mlst}"
        echo -e "\e[32m MLST complete: ${expected_mlst} \e[0m"
    fi

    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m KLEBORATE: ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    expected_kleborate="${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv"
    if [ -s "${expected_kleborate}" ]; then
        echo -e "\e[32m Kleborate output already exists -> skipping \e[0m"
        log_failure "Phase4_Kleborate" "ALL" "SKIPPED_EXISTS" "Output ${expected_kleborate} already exists"
    elif compgen -G "${CHECKM_INPUTS}/*.fasta" > /dev/null; then
        conda run -n BPtyping kleborate \
            --threads "${threads}" \
            --preset kpsc \
            --assemblies "${CHECKM_INPUTS}/"*.fasta \
            --outfile "${expected_kleborate}" 2>/dev/null || true
        echo -e "\e[32m Kleborate complete: ${expected_kleborate} \e[0m"
    fi

    #################################################
    # PHASE 5: ANNOTATION & AMR/VIRULENCE PROFILING #
    #################################################

for sample_id in "${sample_ids[@]}"; do
    sample_name="${sample_id}"
    filtered_Kleb_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_Kleb_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Assembly missing -> skipping annotation \e[0m"
        log_failure "Phase5_Annotation" "${sample_name}" "INPUT_MISSING" "Assembly ${filtered_Kleb_assembly} missing"
        continue
    fi

    ANNOTATION_KLEB_BAKTA="${ANNOTATION_KLEB_BAKTA_BASE}/${sample_name}.bakta"
    ANNOTATION_KLEB_AMRFINDER="${ANNOTATION_KLEB_AMRFINDER_BASE}/${sample_name}.amrfinder"
    EGGNOG_KLEB="${ANNOTATION_KLEB_EGGNOG_BASE}/${sample_name}.eggnog"
    ABRICATE_OUT="${ABRICATE_KLEB}/${sample_name}.abricate"

    mkdir -p "${ANNOTATION_KLEB_BAKTA}"
    mkdir -p "${ANNOTATION_KLEB_AMRFINDER}"
    mkdir -p "${EGGNOG_KLEB}"
    mkdir -p "${ABRICATE_OUT}"

    echo -e "\e[31m =================== \e[0m"
    echo -e "\e[31m BAKTA: ${sample_name} \e[0m"
    echo -e "\e[31m =================== \e[0m"

    bakta_faa="${ANNOTATION_KLEB_BAKTA}/${sample_name}.bakta.faa"
    bakta_gff="${ANNOTATION_KLEB_BAKTA}/${sample_name}.bakta.gff3"

    if [ -s "${bakta_faa}" ] && [ -s "${bakta_gff}" ]; then
        echo -e "\e[32m   [${sample_name}] Bakta annotation already exists -> skipping \e[0m"
        log_failure "Phase5_Bakta" "${sample_name}" "SKIPPED_EXISTS" "Output ${bakta_faa} already exists"
    else
        if ! conda run -n BPannotation bakta \
            --threads "${threads}" \
            --force \
            --db "${bakta_db}" \
            --genus Klebsiella \
            --species pneumoniae \
            --tmp-dir "${ANNOTATION_KLEB_BAKTA}" \
            --output "${ANNOTATION_KLEB_BAKTA}" \
            --prefix "${sample_name}.bakta" \
            "${filtered_Kleb_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Bakta execution failed \e[0m"
            log_failure "Phase5_Bakta" "${sample_name}" "EXECUTION_FAILED" "bakta non-zero exit status"
        else
            echo -e "\e[32m Bakta complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m AMRFINDER: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    expected_amrfinder="${ANNOTATION_KLEB_AMRFINDER}/${sample_name}.amrfinder.tsv"
    if [ -s "${expected_amrfinder}" ]; then
        echo -e "\e[32m   [${sample_name}] AMRFinderPlus output already exists -> skipping \e[0m"
        log_failure "Phase5_AMRFinder" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_amrfinder} already exists"
    elif [ ! -s "${bakta_faa}" ] || [ ! -s "${bakta_gff}" ]; then
        echo -e "\e[31m   [${sample_name}] Bakta FAA or GFF missing -> skipping AMRFinderPlus \e[0m"
        log_failure "Phase5_AMRFinder" "${sample_name}" "INPUT_MISSING" "Bakta FAA/GFF missing for ${sample_name}"
    else
        if ! conda run -n ncbi amrfinder \
            --threads "${threads}" \
            --plus \
            --protein "${bakta_faa}" \
            --gff "${bakta_gff}" \
            --annotation_format bakta \
            --organism Klebsiella_pneumoniae \
            --database "${amrfinder_db}" \
            --mutation_all "${ANNOTATION_KLEB_AMRFINDER}/${sample_name}.mutation_all.amrfinder.tsv" \
            --output "${expected_amrfinder}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: AMRFinder execution failed \e[0m"
            log_failure "Phase5_AMRFinder" "${sample_name}" "EXECUTION_FAILED" "amrfinder non-zero exit status"
        else
            echo -e "\e[32m AMRFinder complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ============================ \e[0m"
    echo -e "\e[31m ABRICATE: ${sample_name} \e[0m"
    echo -e "\e[31m ============================ \e[0m"

    for db in vfdb card ncbi resfinder; do
        expected_abricate="${ABRICATE_OUT}/${sample_name}.abricate.${db}.tsv"
        if [ -s "${expected_abricate}" ]; then
            echo -e "\e[32m   [${sample_name}] ABRICATE ${db} already exists -> skipping \e[0m"
            log_failure "Phase5_ABRICATE_${db}" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_abricate} already exists"
        else
            if ! conda run -n BPannotation abricate \
                --db "${db}" \
                --threads "${threads}" \
                --minid 80 \
                --mincov 80 \
                "${filtered_Kleb_assembly}" \
                > "${expected_abricate}"; then
                echo -e "\e[31m   [${sample_name}] ERROR: ABRICATE ${db} failed \e[0m"
                log_failure "Phase5_ABRICATE_${db}" "${sample_name}" "EXECUTION_FAILED" "abricate --db ${db} failed"
            else
                echo -e "\e[32m ABRICATE ${db} done for ${sample_name} \e[0m"
            fi
        fi
    done

    echo -e "\e[31m =========================== \e[0m"
    echo -e "\e[31m EGGNOG-MAPPER: ${sample_name} \e[0m"
    echo -e "\e[31m =========================== \e[0m"

    expected_eggnog="${EGGNOG_KLEB}/${sample_name}.eggnog.emapper.annotations"
    if [ -s "${expected_eggnog}" ]; then
        echo -e "\e[32m   [${sample_name}] EggNOG-mapper output already exists -> skipping \e[0m"
        log_failure "Phase5_EggNOG" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_eggnog} already exists"
    elif [ ! -s "${bakta_faa}" ]; then
        echo -e "\e[31m   [${sample_name}] Bakta FAA missing -> skipping EggNOG \e[0m"
        log_failure "Phase5_EggNOG" "${sample_name}" "INPUT_MISSING" "Bakta FAA missing for ${sample_name}"
    else
        conda run -n BPannotation emapper.py \
            --cpu "${threads}" \
            --dbmem \
            --data_dir "${eggnog_db}" \
            --temp_dir "${EGGNOG_KLEB}" \
            --output_dir "${EGGNOG_KLEB}" \
            --no_annot \
            --override \
            -i "${bakta_faa}" \
            -o "${sample_name}.eggnog" 2>/dev/null || true

        if [ -f "${EGGNOG_KLEB}/${sample_name}.eggnog.emapper.seed_orthologs" ]; then
            conda run -n BPannotation emapper.py \
                --cpu "${threads}" \
                --data_dir "${eggnog_db}" \
                --temp_dir "${EGGNOG_KLEB}" \
                --output_dir "${EGGNOG_KLEB}" \
                --annotate_hits_table "${EGGNOG_KLEB}/${sample_name}.eggnog.emapper.seed_orthologs" \
                -m no_search \
                --excel \
                --tax_scope 1236 \
                --override \
                -o "${sample_name}.eggnog" 2>/dev/null || true
            echo -e "\e[32m EggNOG-mapper completed for sample ${sample_name} \e[0m"
        fi
    fi
done

    echo -e "\e[31m ============================================== \e[0m"
    echo -e "\e[31m ABRICATE: SUMMARIZE ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ============================================== \e[0m"

    for db in vfdb card ncbi resfinder; do
        conda run -n BPannotation abricate --summary \
            "${ABRICATE_KLEB}/"*".abricate/"*".abricate.${db}.tsv" \
            > "${ABRICATE_KLEB}/all_Kleb.abricate.${db}.summary.tsv" 2>/dev/null || true
        echo -e "\e[32m ABRICATE summary written: all_Kleb.abricate.${db}.summary.tsv \e[0m"
    done

    #################################################
    # PIPELINE EXECUTION SUMMARY                    #
    #################################################

    echo ""
    echo -e "\e[32m ================================================================ \e[0m"
    echo -e "\e[32m K. PNEUMONIAE ASSEMBLY & ANNOTATION (BATCH 2) COMPLETE â€” $(date) \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
    echo ""
    echo -e "\e[32m -- OUTPUT DIRECTORIES ------------------------------------------ \e[0m"
    echo -e "\e[32m  Preprocessing       : ${PREPROCESSING_KLEB} \e[0m"
    echo -e "\e[32m  Assemblies          : ${ASSEMBLY_KLEB} \e[0m"
    echo -e "\e[32m  Annotation          : ${ANNOTATION_KLEB} \e[0m"
    echo -e "\e[32m  Genome completeness : ${ANNOTATION_KLEB_CHECKM} \e[0m"
    echo -e "\e[32m  Sequence typing     : ${ANNOTATION_KLEB_MLST} \e[0m"
    echo -e "\e[32m  Serotype/Virulence  : ${KLEBORATE_KLEB} \e[0m"
    echo -e "\e[32m  AMR (AMRFinderPlus) : ${ANNOTATION_KLEB_AMRFINDER_BASE} \e[0m"
    echo -e "\e[32m  AMR/Virulence (ABR) : ${ABRICATE_KLEB} \e[0m"
    echo ""
    echo -e "\e[32m -- LOG FILES --------------------------------------------------- \e[0m"
    echo -e "\e[32m  Full execution log     : ${LOG} \e[0m"

    fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING' "${FAIL_LOG}" 2>/dev/null || echo 0)
    skip_count=$(grep -c '\[SKIPPED' "${FAIL_LOG}" 2>/dev/null || echo 0)

    if [ "${fail_count}" -gt 0 ]; then
        echo -e "\e[31m  Failure/Issues log     : ${FAIL_LOG} (${fail_count} failures detected!) \e[0m"
        echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which samples failed and why. \e[0m"
    else
        echo -e "\e[32m  Failure/Issues log     : ${FAIL_LOG} (0 errors recorded) \e[0m"
    fi
    echo -e "\e[32m  Skipped checkpoints    : ${skip_count} records \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
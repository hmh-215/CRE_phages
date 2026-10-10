#!/bin/bash
# AMR profiling, virulence, MGE, plasmid assignment, chimera detection, and reassembly of K. pneumoniae isolates (Batch 3 - Hospital 108)
# Workflow:
#   Phase 1 : Extended AMR profiling (RGI CARD contig & protein + ABRICATE multi-database screening)
#   Phase 2 : Mobile Genetic Elements (ISEScan IS elements, IntegronFinder integrons, ICEfinder preparation)
#   Phase 3 : Genome re-annotation (Prokka query isolates & K. pneumoniae NTUH-K2004 reference)
#   Phase 4 : Pathogenicity islands (GIPSy2 vs NTUH-K2004), prophages (PhiSpy), & virulence typing (Kleborate)
#   Phase 5 : Plasmid assignment & reconstruction (Platon, PlasmidFinder, MOB-suite)
#   Phase 6 : Chimera detection & read coverage breakpoints (BWA-MEM, samtools, detect_coverage_breakpoints.py)
#   Phase 7 : Re-assembly & QC (Unicycler, QUAST, CheckM lineage_wf)
# Building No. 2 (Batch 3)
set -euo pipefail

################
# GLOBAL SETUP #
################

# Base directory paths (with fallback handling across naming variants)
BP3_PATH="/storage/student9/projects/bacteria_phages_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage_3" ] && \
    BP3_PATH="/storage/student9/projects/bacteria_phage_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/bacteria_phages_3" ] && \
    BP3_PATH="/storage/student9/bacteria_phages_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/bacteria_phage_3" ] && \
    BP3_PATH="/storage/student9/bacteria_phage_3"

REF_PATH="/storage/student9/references"
SAMPLE_PATH="${BP3_PATH}"
WORK_PATH="${SAMPLE_PATH}/Klebsiella_pneumoniae"
PREPROCESSING_KLEB="${WORK_PATH}/preprocessing_Kleb"
CHECKM_INPUTS="${WORK_PATH}/annotation_Kleb/checkm/inputs_Kleb"
ANNOTATION_KLEB="${WORK_PATH}/annotation_Kleb"

# Tunable execution parameters
threads=16

# Reference genome for comparative analyses (GIPSy2 subtraction)
# K. pneumoniae reference: NTUH-K2004
Kleb_ref="/storage/student9/references/reference_genomes/Klebsiella_pneumoniae/ncbi_dataset/data/GCF_000009885.1/GCF_000009885.1_ASM988v1_genomic.fna"

# Batch 3 sample IDs
sample_ids=(WS2762607A54)
SAMPLE_PREFIX=""

# Paths to external custom tools & scripts
gipsy2="/storage/student9/tools/gipsy/gipsy/gipsy2"
detect_coverage_breakpoints="/storage/student9/tools/detect_coverage_breakpoints.py"

# Paths to reference databases
platon_db="${REF_PATH}/platon_db/db/"
plasmidfinder_db="${REF_PATH}/plasmidfinder_db"
card_db="${REF_PATH}/card_db"

# Per-tool workflow output directories
AMR_KLEB="${WORK_PATH}/amr_Kleb"
RGI_KLEB="${AMR_KLEB}/rgi"
ABRICATE_KLEB="${AMR_KLEB}/abricate"

MGE_KLEB="${WORK_PATH}/mge_Kleb"
ISESCAN_KLEB="${MGE_KLEB}/isescan"
INTEGRON_KLEB="${MGE_KLEB}/integron"
ICEFINDER_KLEB="${MGE_KLEB}/icefinder"

PROKKA_KLEB="${ANNOTATION_KLEB}/prokka"

PAI_KLEB="${WORK_PATH}/pai_Kleb"
GIPSY_KLEB="${PAI_KLEB}/gipsy"
KLEBORATE_KLEB="${PAI_KLEB}/kleborate"
PHISPY_KLEB="${PAI_KLEB}/phispy"

PLASMID_KLEB="${WORK_PATH}/plasmid_Kleb"
PLATON_KLEB="${PLASMID_KLEB}/platon"
PLASMIDFINDER_KLEB="${PLASMID_KLEB}/plasmidfinder"
MOBSUITE_KLEB="${PLASMID_KLEB}/mobsuite"

CHIMERA_KLEB="${WORK_PATH}/chimera_Kleb"
REASSEMBLY_KLEB="${WORK_PATH}/reassembly_Kleb"
REASSEMBLY_QUAST_INPUTS_KLEB="${REASSEMBLY_KLEB}/quast_inputs"

# Create all workflow directories up front
mkdir -p "${WORK_PATH}"
mkdir -p "${RGI_KLEB}"
mkdir -p "${ABRICATE_KLEB}"
mkdir -p "${ISESCAN_KLEB}"
mkdir -p "${INTEGRON_KLEB}"
mkdir -p "${ICEFINDER_KLEB}"
mkdir -p "${PROKKA_KLEB}/prokka_inputs"
mkdir -p "${GIPSY_KLEB}"
mkdir -p "${KLEBORATE_KLEB}"
mkdir -p "${PHISPY_KLEB}"
mkdir -p "${PLATON_KLEB}"
mkdir -p "${PLASMIDFINDER_KLEB}"
mkdir -p "${MOBSUITE_KLEB}"
mkdir -p "${CHIMERA_KLEB}"
mkdir -p "${REASSEMBLY_KLEB}"
mkdir -p "${REASSEMBLY_QUAST_INPUTS_KLEB}"
mkdir -p "${REASSEMBLY_KLEB}/checkm_inputs"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/Kleb_AMRprofiling_reassembly_3.log"
FAIL_LOG="${WORK_PATH}/Kleb_AMRprofiling_reassembly_3.failed_skipped.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================" >> "${FAIL_LOG}"
echo " Failure & Skip Log — Started: $(date)" >> "${FAIL_LOG}"
echo "======================================================" >> "${FAIL_LOG}"

# Helper function to record skipped or failed steps with timestamp
log_failure() {
    local phase="$1"
    local sample="$2"
    local status="$3" # e.g. "SKIPPED_EXISTS", "INPUT_MISSING", "EXECUTION_FAILED", "OUTPUT_MISSING"
    local reason="$4"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${phase}] [${sample}] [${status}] ${reason}" >> "${FAIL_LOG}"
}

echo "================================================================"
echo " K. pneumoniae Extended Profiling & Reassembly (Batch 3) - Started: $(date)"
echo " Total samples in cohort : ${#sample_ids[@]}"
echo " Reference genome        : ${Kleb_ref}"
echo " Full log                : ${LOG}"
echo " Failure & skip log      : ${FAIL_LOG}"
echo "================================================================"

    #################################################
    # PHASE 1: EXTENDED AMR PROFILING (RGI & CARD)  #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 1: EXTENDED AMR PROFILING (RGI & ABRICATE)      \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m RGI: LOAD CARD DATABASE \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    # Load CARD reference ontology once before per-sample execution
    if [ -s "${card_db}/card.json" ]; then
        if ! conda run -n rgi_env rgi load \
            --card_json "${card_db}/card.json" \
            --local; then
            echo -e "\e[31m WARNING: RGI CARD load returned non-zero exit code \e[0m"
            log_failure "Phase1_RGI_load" "ALL" "EXECUTION_FAILED" "Failed to load local CARD JSON"
        else
            echo -e "\e[32m CARD database loaded successfully \e[0m"
        fi
    else
        echo -e "\e[31m WARNING: card.json not found at ${card_db}/card.json \e[0m"
        log_failure "Phase1_RGI_load" "ALL" "INPUT_MISSING" "card.json missing at ${card_db}/card.json"
    fi

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    bakta_faa="${ANNOTATION_KLEB}/bakta/${sample_name}.bakta/${sample_name}.bakta.faa"

    # Checkpoint: Validate input assembly
    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] WARNING: Assembly not found: ${filtered_assembly} - skipping \e[0m"
        log_failure "Phase1_AMR" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    RGI_OUT="${RGI_KLEB}/${sample_name}.rgi"
    ABRICATE_OUT="${ABRICATE_KLEB}/${sample_name}.abricate"
    mkdir -p "${RGI_OUT}"
    mkdir -p "${ABRICATE_OUT}"

    echo -e "\e[31m ================= \e[0m"
    echo -e "\e[31m RGI: ${sample_name} \e[0m"
    echo -e "\e[31m ================= \e[0m"

    # Contig-based RGI: identifies full resistance genes from assembly contigs
    expected_rgi_contig="${RGI_OUT}/${sample_name}.rgi.txt"
    if [ -s "${expected_rgi_contig}" ]; then
        echo -e "\e[32m   [${sample_name}] Contig RGI already completed - skipping \e[0m"
        log_failure "Phase1_RGI_contig" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_rgi_contig} already exists"
    else
        if ! conda run -n rgi_env rgi main \
            --input_sequence "${filtered_assembly}" \
            --output_file "${RGI_OUT}/${sample_name}.rgi" \
            --input_type contig \
            --alignment_tool BLAST \
            --include_loose \
            --include_nudge \
            --num_threads "${threads}" \
            --clean \
            --local; then
            echo -e "\e[31m   [${sample_name}] ERROR: Contig RGI execution failed \e[0m"
            log_failure "Phase1_RGI_contig" "${sample_name}" "EXECUTION_FAILED" "RGI main contig non-zero exit code"
        elif [ ! -s "${expected_rgi_contig}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: Contig RGI expected output missing \e[0m"
            log_failure "Phase1_RGI_contig" "${sample_name}" "OUTPUT_MISSING" "File ${expected_rgi_contig} not found"
        else
            echo -e "\e[32m Contig RGI complete for ${sample_name} -> ${expected_rgi_contig} \e[0m"
        fi
    fi

    # Protein-based RGI: sensitive detection of point mutations (e.g. gyrA, parC)
    expected_rgi_prot="${RGI_OUT}/${sample_name}.rgi.protein.txt"
    if [ -s "${expected_rgi_prot}" ]; then
        echo -e "\e[32m   [${sample_name}] Protein RGI already completed - skipping \e[0m"
        log_failure "Phase1_RGI_protein" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_rgi_prot} already exists"
    elif [ -s "${bakta_faa}" ]; then
        if ! conda run -n rgi_env rgi main \
            --input_sequence "${bakta_faa}" \
            --output_file "${RGI_OUT}/${sample_name}.rgi.protein" \
            --input_type protein \
            --alignment_tool BLAST \
            --include_loose \
            --num_threads "${threads}" \
            --clean \
            --local; then
            echo -e "\e[31m   [${sample_name}] ERROR: Protein RGI execution failed \e[0m"
            log_failure "Phase1_RGI_protein" "${sample_name}" "EXECUTION_FAILED" "RGI main protein non-zero exit code"
        elif [ ! -s "${expected_rgi_prot}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: Protein RGI expected output missing \e[0m"
            log_failure "Phase1_RGI_protein" "${sample_name}" "OUTPUT_MISSING" "File ${expected_rgi_prot} not found"
        else
            echo -e "\e[32m Protein RGI complete for ${sample_name} -> ${expected_rgi_prot} \e[0m"
        fi
    else
        echo -e "\e[32m   [${sample_name}] Bakta protein FAA not found - skipping protein RGI \e[0m"
        log_failure "Phase1_RGI_protein" "${sample_name}" "INPUT_MISSING" "Bakta protein FAA missing: ${bakta_faa}"
    fi

    echo -e "\e[31m ====================== \e[0m"
    echo -e "\e[31m ABRICATE: ${sample_name} \e[0m"
    echo -e "\e[31m ====================== \e[0m"

    # Rapid screening across curated AMR and virulence databases
    for db in card ncbi resfinder argannot vfdb; do
        abricate_tsv="${ABRICATE_OUT}/${sample_name}.abricate.${db}.tsv"
        if [ -s "${abricate_tsv}" ]; then
            echo -e "\e[32m   [${sample_name}] ABRICATE ${db} already completed - skipping \e[0m"
            log_failure "Phase1_ABRI_${db}" "${sample_name}" "SKIPPED_EXISTS" "Output ${abricate_tsv} already exists"
        else
            if ! conda run -n BPannotation abricate \
                --db "${db}" \
                --threads "${threads}" \
                --minid 80 \
                --mincov 80 \
                "${filtered_assembly}" \
                > "${abricate_tsv}"; then
                echo -e "\e[31m   [${sample_name}] ERROR: ABRICATE ${db} execution failed \e[0m"
                log_failure "Phase1_ABRI_${db}" "${sample_name}" "EXECUTION_FAILED" "ABRICATE non-zero exit code"
            elif [ ! -s "${abricate_tsv}" ]; then
                echo -e "\e[31m   [${sample_name}] ERROR: ABRICATE ${db} produced empty output \e[0m"
                log_failure "Phase1_ABRI_${db}" "${sample_name}" "OUTPUT_MISSING" "File ${abricate_tsv} is empty"
            else
                echo -e "\e[32m ABRICATE ${db} complete for ${sample_name} \e[0m"
            fi
        fi
    done
done

    echo -e "\e[31m =========================================== \e[0m"
    echo -e "\e[31m ABRICATE: SUMMARIZE ALL K. PNEUMONIAE SAMPLES\e[0m"
    echo -e "\e[31m =========================================== \e[0m"

    # Aggregate per-database summary tables across all samples
    for db in card ncbi resfinder argannot vfdb; do
        summary_tsv="${ABRICATE_KLEB}/all_Kleb.abricate.${db}.summary.tsv"
        if compgen -G "${ABRICATE_KLEB}/*/*.abricate.${db}.tsv" > /dev/null; then
            conda run -n BPannotation abricate --summary \
                "${ABRICATE_KLEB}/"*".abricate/"*".abricate.${db}.tsv" \
                > "${summary_tsv}"
            echo -e "\e[32m ABRICATE summary written: ${summary_tsv} \e[0m"
        else
            echo -e "\e[33m No ABRICATE ${db} files found to summarize \e[0m"
            log_failure "Phase1_ABRI_Summary" "ALL" "INPUT_MISSING" "No *.abricate.${db}.tsv files found"
        fi
    done

    echo -e "\e[31m ======================================== \e[0m"
    echo -e "\e[31m RGI: HEATMAP OF ALL K. PNEUMONIAE SAMPLES\e[0m"
    echo -e "\e[31m ======================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    json_src="${RGI_KLEB}/${sample_name}.rgi/${sample_name}.rgi.json"
    json_dst="${RGI_KLEB}/${sample_name}rgi.json"
    if [ -f "${json_src}" ] && [ ! -f "${json_dst}" ]; then
        cp -f "${json_src}" "${json_dst}"
    fi
done

    # Generate comparative RGI resistance gene presence/absence heatmap
    expected_heatmap="${RGI_KLEB}/all_Kleb.rgi_heatmap.png"
    if [ -s "${expected_heatmap}" ] || [ -s "${RGI_KLEB}/all_Kleb.rgi_heatmap.eps" ]; then
        echo -e "\e[32m RGI heatmap already exists - skipping \e[0m"
        log_failure "Phase1_RGI_Heatmap" "ALL" "SKIPPED_EXISTS" "Heatmap output already present"
    elif compgen -G "${RGI_KLEB}/*rgi.json" > /dev/null; then
        if ! conda run -n rgi_env rgi heatmap \
            --input "${RGI_KLEB}" \
            --output "${RGI_KLEB}/all_Kleb.rgi_heatmap"; then
            echo -e "\e[31m ERROR: RGI heatmap generation failed \e[0m"
            log_failure "Phase1_RGI_Heatmap" "ALL" "EXECUTION_FAILED" "rgi heatmap non-zero exit code"
        else
            echo -e "\e[32m RGI heatmap generated: ${RGI_KLEB}/all_Kleb.rgi_heatmap \e[0m"
        fi
    else
        echo -e "\e[33m No RGI JSON files found - skipping heatmap \e[0m"
        log_failure "Phase1_RGI_Heatmap" "ALL" "INPUT_MISSING" "No rgi.json files in ${RGI_KLEB}"
    fi

    #############################################
    # PHASE 2: MOBILE GENETIC ELEMENTS (MGEs)   #
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 2: MOBILE GENETIC ELEMENTS (IS, INTEGRONS, ICE) \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_assembly}" ]; then
        log_failure "Phase2_MGE" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    ISESCAN_SAMPLE="${ISESCAN_KLEB}/${sample_name}.isescan"
    INTEGRON_SAMPLE="${INTEGRON_KLEB}/${sample_name}.integron"
    ICEFINDER_SAMPLE="${ICEFINDER_KLEB}/${sample_name}"
    mkdir -p "${ISESCAN_SAMPLE}"
    mkdir -p "${INTEGRON_SAMPLE}"
    mkdir -p "${ICEFINDER_SAMPLE}"

    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m ISESCAN: ${sample_name} \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    # Identify insertion sequences (IS elements, families, terminal repeats)
    if [ -d "${ISESCAN_SAMPLE}/prediction" ] || compgen -G "${ISESCAN_SAMPLE}/*.csv" > /dev/null; then
        echo -e "\e[32m   [${sample_name}] ISEScan prediction already exists - skipping \e[0m"
        log_failure "Phase2_ISEScan" "${sample_name}" "SKIPPED_EXISTS" "ISEScan output already exists"
    else
        if ! conda run -n recombination isescan.py \
            --nthread "${threads}" \
            --seqfile "${filtered_assembly}" \
            --output "${ISESCAN_SAMPLE}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: ISEScan execution failed \e[0m"
            log_failure "Phase2_ISEScan" "${sample_name}" "EXECUTION_FAILED" "isescan.py non-zero exit status"
        else
            echo -e "\e[32m ISEScan complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ============================= \e[0m"
    echo -e "\e[31m INTEGRON FINDER: ${sample_name} \e[0m"
    echo -e "\e[31m ============================= \e[0m"

    # Screen for complete, In0, and CALIN integrons with attC gene cassettes
    INTEGRON_RESULTS="${INTEGRON_SAMPLE}/Results_Integron_Finder_${sample_name}.contigs.filtered"
    if [ -d "${INTEGRON_RESULTS}" ] || compgen -G "${INTEGRON_SAMPLE}/*.integrons" > /dev/null; then
        echo -e "\e[32m   [${sample_name}] IntegronFinder results already exist - skipping \e[0m"
        log_failure "Phase2_Integron" "${sample_name}" "SKIPPED_EXISTS" "IntegronFinder directory already exists"
    else
        if ! conda run -n recombination integron_finder \
            --gbk --pdf \
            --circ \
            --outdir "${INTEGRON_SAMPLE}" \
            "${filtered_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: IntegronFinder execution failed \e[0m"
            log_failure "Phase2_Integron" "${sample_name}" "EXECUTION_FAILED" "integron_finder non-zero exit status"
        else
            # Standardise output naming for downstream visualisation
            if [ -d "${INTEGRON_RESULTS}" ]; then
                for contig_gbk in "${INTEGRON_RESULTS}"/contig_*.gbk; do
                    [ -f "${contig_gbk}" ] || continue
                    c_name=$(basename "${contig_gbk}" .gbk)
                    mv "${contig_gbk}" "${INTEGRON_RESULTS}/${sample_name}.${c_name}.integron.gbk"
                done

                pdf_count=0
                for pdf_file in "${INTEGRON_RESULTS}"/contig_*_*.pdf; do
                    [ -f "${pdf_file}" ] || continue
                    integron_num=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f1 | rev)
                    contig_id=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f2- | rev | sed 's/^contig_//')
                    mv "${pdf_file}" "${INTEGRON_RESULTS}/${sample_name}.contig${contig_id}.integron_${integron_num}.pdf"
                    pdf_count=$(( pdf_count + 1 ))
                done
                echo -e "\e[32m IntegronFinder PDFs structured: ${pdf_count} for ${sample_name} \e[0m"
            fi
            echo -e "\e[32m IntegronFinder complete for ${sample_name} \e[0m"
        fi
    fi

    # Prepare ICEfinder query package (web server submission)
    cp -f "${filtered_assembly}" "${ICEFINDER_SAMPLE}/${sample_name}_for_ICEfinder.fasta"
done

    echo -e "\e[32m ICEfinder queries formatted in: ${ICEFINDER_KLEB}/<sample>/ \e[0m"
    echo -e "\e[32m (Upload at: https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html) \e[0m"

    #############################################
    # PHASE 3: GENOME RE-ANNOTATION (PROKKA)    #
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 3: GENOME RE-ANNOTATION (PROKKA)                \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    renamed_assembly="${PROKKA_KLEB}/prokka_inputs/${sample_name}.renamed.fasta"
    sample_prokka_dir="${PROKKA_KLEB}/${sample_name}.prokka"
    expected_prokka_gbk="${sample_prokka_dir}/${sample_name}.prokka.gbk"

    if [ ! -s "${filtered_assembly}" ]; then
        log_failure "Phase3_Prokka" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    # Clean and standardize contig headers (required for GIPSy & PhiSpy compatibility)
    if [ ! -s "${renamed_assembly}" ]; then
        awk '/^>/ {
            counter++
            print ">contig_" counter
            next
        }
        { print }' "${filtered_assembly}" > "${renamed_assembly}"
        echo -e "\e[32m Renamed headers prepared for ${sample_name} \e[0m"
    fi

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PROKKA: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    if [ -s "${expected_prokka_gbk}" ]; then
        echo -e "\e[32m   [${sample_name}] Prokka annotation already completed - skipping \e[0m"
        log_failure "Phase3_Prokka" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_prokka_gbk} already exists"
    else
        if ! conda run -n BPannotation prokka \
            --force \
            --cpus "${threads}" \
            --genus Klebsiella \
            --species pneumoniae \
            --prefix "${sample_name}.prokka" \
            --outdir "${sample_prokka_dir}" \
            "${renamed_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Prokka annotation failed \e[0m"
            log_failure "Phase3_Prokka" "${sample_name}" "EXECUTION_FAILED" "prokka non-zero exit code"
        elif [ ! -s "${expected_prokka_gbk}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: Prokka expected GBK output missing \e[0m"
            log_failure "Phase3_Prokka" "${sample_name}" "OUTPUT_MISSING" "File ${expected_prokka_gbk} was not created"
        else
            echo -e "\e[32m Prokka annotation complete for ${sample_name} -> ${expected_prokka_gbk} \e[0m"
        fi
    fi
done

    # Annotate K. pneumoniae reference NTUH-K2004 if not already annotated
    REF_PROKKA_DIR="${PROKKA_KLEB}/NTUH-K2004.prokka"
    expected_ref_gbk="${REF_PROKKA_DIR}/NTUH-K2004.prokka.gbk"
    echo -e "\e[31m ============================== \e[0m"
    echo -e "\e[31m PROKKA: REF NTUH-K2004         \e[0m"
    echo -e "\e[31m ============================== \e[0m"

    if [ -s "${Kleb_ref}" ]; then
        if [ -s "${expected_ref_gbk}" ]; then
            echo -e "\e[32m Reference NTUH-K2004 Prokka GBK already exists - skipping \e[0m"
            log_failure "Phase3_Prokka_Ref" "NTUH-K2004" "SKIPPED_EXISTS" "Output ${expected_ref_gbk} already exists"
        else
            if ! conda run -n BPannotation prokka \
                --force \
                --cpus "${threads}" \
                --genus Klebsiella \
                --species pneumoniae \
                --prefix "NTUH-K2004.prokka" \
                --outdir "${REF_PROKKA_DIR}" \
                "${Kleb_ref}"; then
                echo -e "\e[31m ERROR: Reference Prokka annotation failed \e[0m"
                log_failure "Phase3_Prokka_Ref" "NTUH-K2004" "EXECUTION_FAILED" "prokka non-zero exit on reference"
            else
                echo -e "\e[32m Reference NTUH-K2004 Prokka annotation complete -> ${expected_ref_gbk} \e[0m"
            fi
        fi
    else
        echo -e "\e[31m WARNING: K. pneumoniae reference genome not found at: ${Kleb_ref} \e[0m"
        log_failure "Phase3_Prokka_Ref" "NTUH-K2004" "INPUT_MISSING" "Reference FASTA not found at ${Kleb_ref}"
    fi

    #################################################
    # PHASE 4: PATHOGENICITY ISLANDS & VIRULENCE    #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 4: PATHOGENICITY ISLANDS (GIPSY2), PHISPY & KLEB\e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    # Export LD_LIBRARY_PATH for gipsy2 binary library dependencies
    export LD_LIBRARY_PATH="/storage/student9/miniconda3/envs/gipsy_env/lib:${LD_LIBRARY_PATH:-}"
    Kleb_ref_gbk="${REF_PROKKA_DIR}/NTUH-K2004.prokka.gbk"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    prokka_gbk="${PROKKA_KLEB}/${sample_name}.prokka/${sample_name}.prokka.gbk"
    gipsy_out="${GIPSY_KLEB}/${sample_name}.gipsy2"
    phispy_out_dir="${PHISPY_KLEB}/${sample_name}.phispy"

    if [ ! -s "${prokka_gbk}" ]; then
        echo -e "\e[31m   [${sample_name}] Prokka GBK missing - skipping PAI/prophage analysis \e[0m"
        log_failure "Phase4_PAI" "${sample_name}" "INPUT_MISSING" "Prokka GBK ${prokka_gbk} missing"
        continue
    fi

    mkdir -p "${gipsy_out}"
    mkdir -p "${phispy_out_dir}"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m GIPSY2: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    # Comparative island prediction subtracting NTUH-K2004 backbone
    expected_gipsy_txt="${gipsy_out}/pathogenicity_islands.txt"
    if [ -s "${Kleb_ref_gbk}" ]; then
        if [ -s "${expected_gipsy_txt}" ]; then
            echo -e "\e[32m   [${sample_name}] GIPSy2 output already exists - skipping \e[0m"
            log_failure "Phase4_GIPSy2" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_gipsy_txt} already exists"
        else
            if ! conda run -n gipsy_env "${gipsy2}" \
                -q "${prokka_gbk}" \
                -s "${Kleb_ref_gbk}" \
                -o "${gipsy_out}" \
                -res -vir -met \
                -k fisher \
                --force; then
                echo -e "\e[31m   [${sample_name}] ERROR: GIPSy2 execution failed \e[0m"
                log_failure "Phase4_GIPSy2" "${sample_name}" "EXECUTION_FAILED" "gipsy2 non-zero exit status"
            else
                echo -e "\e[32m GIPSy2 comparative profiling completed for ${sample_name} \e[0m"
            fi
        fi
    else
        echo -e "\e[31m   [${sample_name}] Reference GBK ${Kleb_ref_gbk} not available - skipping GIPSy2 \e[0m"
        log_failure "Phase4_GIPSy2" "${sample_name}" "INPUT_MISSING" "Reference GBK ${Kleb_ref_gbk} missing"
    fi

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PHISPY: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    # Identify prophages integrated in the genome
    expected_phispy_tsv="${phispy_out_dir}/prophage.tsv"
    if [ -s "${expected_phispy_tsv}" ]; then
        echo -e "\e[32m   [${sample_name}] PhiSpy prediction already exists - skipping \e[0m"
        log_failure "Phase4_PhiSpy" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_phispy_tsv} already exists"
    else
        if ! conda run -n recombination PhiSpy.py \
            "${prokka_gbk}" \
            -o "${phispy_out_dir}" \
            --output_choice 7 \
            --threads "${threads}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: PhiSpy execution failed \e[0m"
            log_failure "Phase4_PhiSpy" "${sample_name}" "EXECUTION_FAILED" "PhiSpy.py non-zero exit status"
        else
            echo -e "\e[32m PhiSpy complete for ${sample_name} -> ${expected_phispy_tsv} \e[0m"
        fi
    fi
done

    # Merge prophage coordinates and records across all cohort isolates
    echo -e "\e[31m ============================ \e[0m"
    echo -e "\e[31m MERGING PHISPY: ALL SAMPLES \e[0m"
    echo -e "\e[31m ============================ \e[0m"

    MERGED_PHISPY="${PHISPY_KLEB}/all_Kleb.phispy.merged.tsv"
    : > "${MERGED_PHISPY}"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    phispy_tsv="${PHISPY_KLEB}/${sample_name}.phispy/prophage.tsv"

    if [ -s "${phispy_tsv}" ]; then
        awk -v s="${sample_name}" 'BEGIN{FS=OFS="\t"} {print s, $0}' "${phispy_tsv}" >> "${MERGED_PHISPY}"
        echo -e "\e[32m   ${sample_name}: $(wc -l < "${phispy_tsv}") prophage(s) merged \e[0m"
    fi
done

    echo -e "\e[32m Merged PhiSpy table: ${MERGED_PHISPY} (Total: $(wc -l < "${MERGED_PHISPY}") records) \e[0m"

    # Kleborate: dedicated Klebsiella virulence, capsular serotyping and resistance profiling
    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m KLEBORATE: ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    if [ -s "${KLEBORATE_KLEB}/klebsiella_pneumo_complex_output.txt" ] || [ -s "${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv" ]; then
        echo -e "\e[32m Kleborate output already exists - skipping \e[0m"
        log_failure "Phase4_Kleborate" "ALL_COHORT" "SKIPPED_EXISTS" "Kleborate output table already exists"
    elif compgen -G "${CHECKM_INPUTS}/*.fasta" > /dev/null; then
        if ! conda run -n BPtyping kleborate \
            --threads "${threads}" \
            --preset kpsc \
            --assemblies "${CHECKM_INPUTS}/"*.fasta \
            --outdir "${KLEBORATE_KLEB}"; then
            echo -e "\e[31m ERROR: Kleborate execution failed \e[0m"
            log_failure "Phase4_Kleborate" "ALL_COHORT" "EXECUTION_FAILED" "kleborate command exited non-zero"
        else
            if [ -s "${KLEBORATE_KLEB}/klebsiella_pneumo_complex_output.txt" ]; then
                cp -f "${KLEBORATE_KLEB}/klebsiella_pneumo_complex_output.txt" "${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv"
            fi
            echo -e "\e[32m Kleborate complete -> ${KLEBORATE_KLEB}/ \e[0m"
        fi
    else
        echo -e "\e[31m ERROR: No assembly fasta files in ${CHECKM_INPUTS} \e[0m"
        log_failure "Phase4_Kleborate" "ALL_COHORT" "INPUT_MISSING" "No fasta files in ${CHECKM_INPUTS}"
    fi

    #################################################
    # PHASE 5: PLASMID IDENTIFICATION & TYPING      #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 5: PLASMID ASSIGNMENT (PLATON, MOB-SUITE, INC)  \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_assembly}" ]; then
        log_failure "Phase5_Plasmid" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    PLATON_OUT="${PLASMID_KLEB}/${sample_name}.platon"
    PLASMIDFINDER_OUT="${PLASMID_KLEB}/${sample_name}.plasmidfinder"
    MOBSUITE_OUT="${PLASMID_KLEB}/${sample_name}.mobsuite"

    mkdir -p "${PLATON_OUT}"
    mkdir -p "${PLASMIDFINDER_OUT}"
    mkdir -p "${MOBSUITE_OUT}"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PLATON: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    # Classify contigs as chromosomal vs plasmid based on marker HMMs and replication proteins
    expected_platon="${PLATON_OUT}/${sample_name}.chromosome.fasta"
    if [ -s "${expected_platon}" ]; then
        echo -e "\e[32m   [${sample_name}] Platon classification already exists - skipping \e[0m"
        log_failure "Phase5_Platon" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_platon} already exists"
    else
        if ! conda run -n plasmid platon \
            --db "${platon_db}" \
            --output "${PLATON_OUT}" \
            --prefix "${sample_name}" \
            --mode sensitivity \
            --threads "${threads}" \
            "${filtered_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Platon execution failed \e[0m"
            log_failure "Phase5_Platon" "${sample_name}" "EXECUTION_FAILED" "platon non-zero exit status"
        else
            echo -e "\e[32m Platon complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m =========================== \e[0m"
    echo -e "\e[31m PLASMIDFINDER: ${sample_name} \e[0m"
    echo -e "\e[31m =========================== \e[0m"

    # Identify incompatibility (Inc) groups and replicon types
    expected_plasmidfinder="${PLASMIDFINDER_OUT}/results_tab.tsv"
    if [ -s "${expected_plasmidfinder}" ]; then
        echo -e "\e[32m   [${sample_name}] PlasmidFinder output already exists - skipping \e[0m"
        log_failure "Phase5_PlasmidFinder" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_plasmidfinder} already exists"
    else
        if ! conda run -n plasmid plasmidfinder.py \
            -i "${filtered_assembly}" \
            -o "${PLASMIDFINDER_OUT}" \
            -p "${plasmidfinder_db}" \
            -l 0.60 \
            -t 0.80 \
            -x; then
            echo -e "\e[31m   [${sample_name}] ERROR: PlasmidFinder execution failed \e[0m"
            log_failure "Phase5_PlasmidFinder" "${sample_name}" "EXECUTION_FAILED" "plasmidfinder.py non-zero exit status"
        else
            echo -e "\e[32m PlasmidFinder complete for ${sample_name} -> ${expected_plasmidfinder} \e[0m"
        fi
    fi

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m MOB-SUITE: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    # Cluster contigs into complete plasmids, identify relaxase and mate-pair formation (mpf)
    expected_mob="${MOBSUITE_OUT}/mobtyper_results.txt"
    if [ -s "${expected_mob}" ]; then
        echo -e "\e[32m   [${sample_name}] MOB-suite results already exist - skipping \e[0m"
        log_failure "Phase5_MOBsuite" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_mob} already exists"
    else
        if ! conda run -n plasmid mob_recon \
            --infile "${filtered_assembly}" \
            --outdir "${MOBSUITE_OUT}" \
            --num_threads "${threads}" \
            --force; then
            echo -e "\e[31m   [${sample_name}] ERROR: MOB-suite execution failed \e[0m"
            log_failure "Phase5_MOBsuite" "${sample_name}" "EXECUTION_FAILED" "mob_recon non-zero exit status"
        else
            echo -e "\e[32m MOB-suite complete for ${sample_name} -> ${expected_mob} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 6: CHIMERA DETECTION & READ COVERAGE    #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 6: READ MAPPING & CHIMERIC BREAKPOINT DETECTION \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    read1t="${PREPROCESSING_KLEB}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_KLEB}/${sample_name}.R2.paired.fastq.gz"

    if [ ! -s "${filtered_assembly}" ] || [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        echo -e "\e[31m   [${sample_name}] Fastq or assembly missing - skipping chimera check \e[0m"
        log_failure "Phase6_Chimera" "${sample_name}" "INPUT_MISSING" "Reads or assembly missing for ${sample_name}"
        continue
    fi

    CHIMERA_OUT="${CHIMERA_KLEB}/${sample_name}.chimera"
    mkdir -p "${CHIMERA_OUT}"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m BWA INDEX: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    if [ ! -s "${filtered_assembly}.bwt" ]; then
        conda run -n mapping bwa index "${filtered_assembly}"
    fi

    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m BWA MEM: ${sample_name} \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    # Map reads, sort and index BAM
    mapped_bam="${CHIMERA_OUT}/${sample_name}.mapped.bam"
    if [ -s "${mapped_bam}" ]; then
        echo -e "\e[32m   [${sample_name}] BAM already mapped - skipping BWA \e[0m"
        log_failure "Phase6_BWA" "${sample_name}" "SKIPPED_EXISTS" "Output ${mapped_bam} already exists"
    else
        conda run -n mapping bwa mem \
            -t "${threads}" \
            "${filtered_assembly}" \
            "${read1t}" "${read2t}" \
            -o "${CHIMERA_OUT}/${sample_name}.tmp.sam" \
            2>> "${CHIMERA_OUT}/${sample_name}.bwa.log"

        conda run -n mapping samtools sort \
            -@ "${threads}" \
            -o "${mapped_bam}" \
            "${CHIMERA_OUT}/${sample_name}.tmp.sam"

        rm -f "${CHIMERA_OUT}/${sample_name}.tmp.sam"
        conda run -n mapping samtools index "${mapped_bam}"
    fi

    # Per-contig summary & per-base coverage profile
    coverage_summary="${CHIMERA_OUT}/${sample_name}.coverage_summary.tsv"
    depth_per_base="${CHIMERA_OUT}/${sample_name}.depth_per_base.tsv"
    flagstat_txt="${CHIMERA_OUT}/${sample_name}.flagstat.txt"

    [ -s "${coverage_summary}" ] || conda run -n mapping samtools coverage "${mapped_bam}" > "${coverage_summary}"
    [ -s "${depth_per_base}" ] || conda run -n mapping samtools depth -a "${mapped_bam}" > "${depth_per_base}"
    [ -s "${flagstat_txt}" ] || conda run -n mapping samtools flagstat "${mapped_bam}" > "${flagstat_txt}"

    echo -e "\e[32m Coverage statistics extracted for ${sample_name} \e[0m"

    # Detect coverage breakpoints (misassembly or chimeric joins)
    if [ -f "${detect_coverage_breakpoints}" ]; then
        breakpoints_report="${CHIMERA_OUT}/${sample_name}.breakpoints_report.tsv"
        if [ -s "${breakpoints_report}" ]; then
            echo -e "\e[32m   [${sample_name}] Breakpoints report already exists - skipping \e[0m"
            log_failure "Phase6_Breakpoints" "${sample_name}" "SKIPPED_EXISTS" "Output ${breakpoints_report} already exists"
        else
            echo -e "\e[32m Running breakpoint detection for ${sample_name}... \e[0m"
            if ! conda run -n ncbi python3 "${detect_coverage_breakpoints}" \
                --input "${depth_per_base}" \
                --coverage-summary "${coverage_summary}" \
                --window 500 \
                --cov-cutoff 0.2 --drop-size 50 \
                --jump-cutoff 3.0 --jump-size 50 \
                --output "${breakpoints_report}"; then
                echo -e "\e[31m   [${sample_name}] ERROR: Breakpoint detection script failed \e[0m"
                log_failure "Phase6_Breakpoints" "${sample_name}" "EXECUTION_FAILED" "detect_coverage_breakpoints.py error"
            else
                echo -e "\e[32m Breakpoints report generated for ${sample_name} -> ${breakpoints_report} \e[0m"
            fi
        fi
    fi
done

    # Summarise mapping rates across cohort
    echo -e "\e[32m Mapping rate summary (all samples): \e[0m"
    for sample_id in "${sample_ids[@]}"; do
        sample_name="${SAMPLE_PREFIX}${sample_id}"
        f_stat="${CHIMERA_KLEB}/${sample_name}.chimera/${sample_name}.flagstat.txt"
        if [ -s "${f_stat}" ]; then
            rate=$(grep "mapped (" "${f_stat}" | head -1 | awk '{print $5}' | tr -d '()')
            echo "  ${sample_name}: ${rate}"
        fi
    done

    #################################################
    # PHASE 7: RE-ASSEMBLY (UNICYCLER) & QC         #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 7: RE-ASSEMBLY (UNICYCLER) & QUALITY QC (CHECKM)\e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    read1t="${PREPROCESSING_KLEB}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_KLEB}/${sample_name}.R2.paired.fastq.gz"

    if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        log_failure "Phase7_Unicycler" "${sample_name}" "INPUT_MISSING" "Paired reads missing for ${sample_name}"
        continue
    fi

    UNICYCLER_OUT="${REASSEMBLY_KLEB}/${sample_name}.unicycler"
    expected_assembly="${UNICYCLER_OUT}/assembly.fasta"
    mkdir -p "${UNICYCLER_OUT}"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m UNICYCLER: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    if [ -s "${expected_assembly}" ]; then
        echo -e "\e[32m   [${sample_name}] Unicycler assembly already exists - skipping \e[0m"
        log_failure "Phase7_Unicycler" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_assembly} already exists"
    else
        if ! conda run -n assembly unicycler \
            -1 "${read1t}" \
            -2 "${read2t}" \
            --mode conservative \
            --threads "${threads}" \
            --out "${UNICYCLER_OUT}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Unicycler execution failed \e[0m"
            log_failure "Phase7_Unicycler" "${sample_name}" "EXECUTION_FAILED" "unicycler exited non-zero"
            continue
        elif [ ! -s "${expected_assembly}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: Unicycler expected assembly.fasta missing \e[0m"
            log_failure "Phase7_Unicycler" "${sample_name}" "OUTPUT_MISSING" "File ${expected_assembly} not found"
            continue
        else
            echo -e "\e[32m Unicycler complete for ${sample_name} -> ${expected_assembly} \e[0m"
        fi
    fi

    # QUAST single-isolate assembly metrics
    if [ -s "${expected_assembly}" ]; then
        [ -d "${UNICYCLER_OUT}/quast" ] || conda run -n assembly quast.py \
            "${expected_assembly}" \
            --threads "${threads}" \
            --output-dir "${UNICYCLER_OUT}/quast"

        # Copy to multi-sample inputs
        cp -f "${expected_assembly}" "${REASSEMBLY_QUAST_INPUTS_KLEB}/${sample_name}.unicycler.fasta"
        cp -f "${expected_assembly}" "${REASSEMBLY_KLEB}/checkm_inputs/${sample_name}.unicycler.fasta"
    fi
done

    # Multi-sample comparative QUAST summary
    echo -e "\e[31m ========================================= \e[0m"
    echo -e "\e[31m QUAST: ALL UNICYCLER ASSEMBLIES (SUMMARY) \e[0m"
    echo -e "\e[31m ========================================= \e[0m"

    if compgen -G "${REASSEMBLY_QUAST_INPUTS_KLEB}/*.fasta" > /dev/null; then
        conda run -n assembly quast.py \
            --threads "${threads}" \
            --output-dir "${REASSEMBLY_KLEB}/all_Kleb.unicycler.quast" \
            "${REASSEMBLY_QUAST_INPUTS_KLEB}/"*.fasta
        echo -e "\e[32m Multi-sample QUAST complete: ${REASSEMBLY_KLEB}/all_Kleb.unicycler.quast \e[0m"
    fi

    # Multi-sample CheckM lineage completeness QC
    echo -e "\e[31m ================================================== \e[0m"
    echo -e "\e[31m CHECKM: ALL UNICYCLER ASSEMBLIES (COMPLETENESS QC)\e[0m"
    echo -e "\e[31m ================================================== \e[0m"

    if compgen -G "${REASSEMBLY_KLEB}/checkm_inputs/*.fasta" > /dev/null; then
        conda run -n BPstructure checkm lineage_wf \
            -t "${threads}" \
            --reduced_tree \
            --pplacer_threads 1 \
            -x fasta \
            --tab_table \
            "${REASSEMBLY_KLEB}/checkm_inputs" \
            "${REASSEMBLY_KLEB}/checkm"

        if [ -f "${REASSEMBLY_KLEB}/checkm/lineage.ms" ]; then
            conda run -n BPstructure checkm qa \
                "${REASSEMBLY_KLEB}/checkm/lineage.ms" \
                "${REASSEMBLY_KLEB}/checkm" \
                -o 2 \
                --tab_table \
                -f "${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv"
            echo -e "\e[32m CheckM QA complete: ${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv \e[0m"
        fi
    fi

    #################################################
    # PIPELINE EXECUTION SUMMARY                    #
    #################################################

    echo ""
    echo -e "\e[32m ================================================================ \e[0m"
    echo -e "\e[32m K. PNEUMONIAE PROFILING & REASSEMBLY (BATCH 3) COMPLETE - $(date)\e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
    echo ""
    echo -e "\e[32m -- OUTPUT DIRECTORIES ------------------------------------------ \e[0m"
    echo -e "\e[32m  AMR (RGI CARD)         : ${RGI_KLEB}/ \e[0m"
    echo -e "\e[32m  AMR (ABRICATE multi-db): ${ABRICATE_KLEB}/ \e[0m"
    echo -e "\e[32m  IS Elements (ISEScan)  : ${ISESCAN_KLEB}/ \e[0m"
    echo -e "\e[32m  Integrons (IntegronF.) : ${INTEGRON_KLEB}/ \e[0m"
    echo -e "\e[32m  ICEfinder FASTA inputs : ${ICEFINDER_KLEB}/ \e[0m"
    echo -e "\e[32m  Prokka Re-annotations  : ${PROKKA_KLEB}/ \e[0m"
    echo -e "\e[32m  PAIs (GIPSy2)          : ${GIPSY_KLEB}/ \e[0m"
    echo -e "\e[32m  Prophages (PhiSpy)     : ${PHISPY_KLEB}/all_Kleb.phispy.merged.tsv \e[0m"
    echo -e "\e[32m  Virulence (Kleborate)  : ${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv \e[0m"
    echo -e "\e[32m  Platon Plasmids        : ${PLATON_KLEB}/ \e[0m"
    echo -e "\e[32m  PlasmidFinder (Inc)    : ${PLASMIDFINDER_KLEB}/ \e[0m"
    echo -e "\e[32m  MOB-suite Plasmids     : ${MOBSUITE_KLEB}/ \e[0m"
    echo -e "\e[32m  Chimera Breakpoint QC  : ${CHIMERA_KLEB}/ \e[0m"
    echo -e "\e[32m  Unicycler Reassembly   : ${REASSEMBLY_KLEB}/ \e[0m"
    echo -e "\e[32m  CheckM QC summary      : ${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv \e[0m"
    echo ""
    echo -e "\e[32m -- LOG FILES --------------------------------------------------- \e[0m"
    echo -e "\e[32m  Full execution log     : ${LOG} \e[0m"

    # Report count of failed or skipped samples recorded in FAIL_LOG
    fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING\]\|\[INPUT_MISSING\]' "${FAIL_LOG}" 2>/dev/null || echo 0)
    skip_count=$(grep -c '\[SKIPPED' "${FAIL_LOG}" 2>/dev/null || echo 0)

    if [ "${fail_count}" -gt 0 ]; then
        echo -e "\e[31m  Failure/Issues log     : ${FAIL_LOG} (${fail_count} failures detected!) \e[0m"
        echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which samples failed and why. \e[0m"
    else
        echo -e "\e[32m  Failure/Issues log     : ${FAIL_LOG} (0 errors recorded) \e[0m"
    fi
    echo -e "\e[32m  Skipped checkpoints    : ${skip_count} records \e[0m"
    echo ""
    echo -e "\e[32m -- SUGGESTED MANUAL & UPSTREAM STEPS --------------------------- \e[0m"
    echo -e "\e[32m  1. Submit formatted ICEfinder FASTA files to: \e[0m"
    echo -e "\e[32m     https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html \e[0m"
    echo -e "\e[32m  2. Open assembly.gfa in Bandage to confirm circular topology \e[0m"
    echo -e "\e[32m  3. Review breakpoints_report.tsv to inspect any chimeric joins \e[0m"
    echo -e "\e[32m  4. Cross-reference Platon + MOB-suite + PlasmidFinder to assign \e[0m"
    echo -e "\e[32m     resistance genes to specific plasmid replicons or chromosomes \e[0m"
    echo -e "\e[32m  5. Check Kleborate output for yersiniabactin/colibactin/aerobactin \e[0m"
    echo -e "\e[32m     virulence loci alongside K/O capsular serotypes \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"

#!/bin/bash
# Escherichia coli AMR profiling, serotyping, virulence, MGE, PAI, and reassembly
# Adapted from Kleb_AMRprofiling_reassembly.bash for E. coli isolates
# Workflow:
#   Phase 1 : Extended AMR profiling (RGI CARD + ABRICATE multi-database screening)
#   Phase 2 : Mobile Genetic Elements (ISEScan IS elements, IntegronFinder integrons, ICEfinder prep)
#   Phase 3 : Genome re-annotation (Prokka query isolates & E. coli K-12 MG1655 reference)
#   Phase 4 : Pathogenicity islands (GIPSy2 vs K-12 MG1655) & prophages (PhiSpy)
#   Phase 5 : E. coli serotyping (ECTyper), phylogrouping (Clermontyping), MLST & virulence (VFDB/AMRFinderPlus)
#   Phase 6 : Plasmid assignment & reconstruction (Platon, PlasmidFinder, MOB-suite)
#   Phase 7 : Chimera detection & read coverage breakpoints (BWA-MEM, samtools, detect_coverage_breakpoints.py)
#   Phase 8 : Re-assembly & QC (Unicycler, QUAST, CheckM lineage_wf)
# Building no. 2 (Batch 2)
set -euo pipefail

################
# GLOBAL SETUP #
################

# Base directory paths
REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages_2"
WORK_PATH="${SAMPLE_PATH}/Escherichia_coli"

PREPROCESSING_ECOLI="${WORK_PATH}/preprocessing_Ecoli"
CHECKM_INPUTS="${WORK_PATH}/annotation_Ecoli/checkm/inputs_Ecoli"
ANNOTATION_ECOLI="${WORK_PATH}/annotation_Ecoli"

# Tunable execution parameters
threads=16

# Reference genome for E. coli comparative analyses (GIPSy2 subtraction)
# Using standard non-pathogenic reference: Escherichia coli str. K-12 substr. MG1655
Ecoli_ref="/storage/student9/references/reference_genomes/Escherichia_coli/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"

#E. coli samples from the ktest submission sheet (Sample # column)
sample_ids=(132895-LJF30331 134194-LJF30332 144342-LJF30333 156589-LJF30334 194083-LJF30335
216381-LJF30336 216607-LJF30337 243686-LJF30338 252837-LJF30339 257625-LJF30340
319707-LJF30341 332413-LJF30342 346152-LJF30343 125919-LJF30344 131695-LJF30345
196951-LJF30346 222949-LJF30347 241208-LJF30348 251713-LJF30349 247233-LJF30350
256139-LJF30351 258569-LJF30352 EM22180-LJF30353 EM10288-LJF30354 EM21266-LJF30360
EM16897-LJF30361-W1 EM10599-LJF30364 EM10610-LJF30365 EM10311-LJF30362)

# Optional sample prefix (set e.g. "EC_" or "WS2762512A" to prepend to sample_ids, or leave empty "" if sample_ids has full names)
SAMPLE_PREFIX=""

# Paths to external custom tools & scripts
gipsy2="/storage/student9/tools/gipsy/gipsy/gipsy2"
detect_coverage_breakpoints="/storage/student9/tools/detect_coverage_breakpoints.py"

# Paths to reference databases
platon_db="${REF_PATH}/platon_db/db/"
plasmidfinder_db="${REF_PATH}/plasmidfinder_db"
card_db="${REF_PATH}/card_db"

# Per-tool workflow output directories
AMR_ECOLI="${WORK_PATH}/amr_Ecoli"
RGI_ECOLI="${AMR_ECOLI}/rgi"
ABRICATE_ECOLI="${AMR_ECOLI}/abricate"

MGE_ECOLI="${WORK_PATH}/mge_Ecoli"
ISESCAN_ECOLI="${MGE_ECOLI}/isescan"
INTEGRON_ECOLI="${MGE_ECOLI}/integron"
ICEFINDER_ECOLI="${MGE_ECOLI}/icefinder"

PROKKA_ECOLI="${ANNOTATION_ECOLI}/prokka"

PAI_ECOLI="${WORK_PATH}/pai_Ecoli"
GIPSY_ECOLI="${PAI_ECOLI}/gipsy"
PHISPY_ECOLI="${PAI_ECOLI}/phispy"

TYPING_ECOLI="${WORK_PATH}/typing_Ecoli"
ECTYPER_ECOLI="${TYPING_ECOLI}/ectyper"
CLERMONT_ECOLI="${TYPING_ECOLI}/clermontyping"
MLST_ECOLI="${TYPING_ECOLI}/mlst"
VIRULENCE_ECOLI="${TYPING_ECOLI}/virulence"

PLASMID_ECOLI="${WORK_PATH}/plasmid_Ecoli"
PLATON_ECOLI="${PLASMID_ECOLI}/platon"
PLASMIDFINDER_ECOLI="${PLASMID_ECOLI}/plasmidfinder"
MOBSUITE_ECOLI="${PLASMID_ECOLI}/mobsuite"

CHIMERA_ECOLI="${WORK_PATH}/chimera_Ecoli"
REASSEMBLY_ECOLI="${WORK_PATH}/reassembly_Ecoli"
REASSEMBLY_QUAST_INPUTS_ECOLI="${REASSEMBLY_ECOLI}/quast_inputs"

# Create all workflow directories up front
mkdir -p "${RGI_ECOLI}"
mkdir -p "${ABRICATE_ECOLI}"
mkdir -p "${ISESCAN_ECOLI}"
mkdir -p "${INTEGRON_ECOLI}"
mkdir -p "${ICEFINDER_ECOLI}"
mkdir -p "${PROKKA_ECOLI}/prokka_inputs"
mkdir -p "${GIPSY_ECOLI}"
mkdir -p "${PHISPY_ECOLI}"
mkdir -p "${ECTYPER_ECOLI}"
mkdir -p "${CLERMONT_ECOLI}"
mkdir -p "${MLST_ECOLI}"
mkdir -p "${VIRULENCE_ECOLI}"
mkdir -p "${PLATON_ECOLI}"
mkdir -p "${PLASMIDFINDER_ECOLI}"
mkdir -p "${MOBSUITE_ECOLI}"
mkdir -p "${CHIMERA_ECOLI}"
mkdir -p "${REASSEMBLY_ECOLI}"
mkdir -p "${REASSEMBLY_QUAST_INPUTS_ECOLI}"
mkdir -p "${REASSEMBLY_ECOLI}/checkm_inputs"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/Ecoli_AMRprofiling_reassembly_2.log"
FAIL_LOG="${WORK_PATH}/Ecoli_AMRprofiling_reassembly_2.failed_skipped.log"
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
echo " E. coli Extended Profiling & Reassembly Pipeline - Started: $(date)"
echo " Total samples in cohort : ${#sample_ids[@]}"
echo " Reference genome        : ${Ecoli_ref}"
echo " Full log                : ${LOG}"
echo " Failure & skip log      : ${FAIL_LOG}"
echo "================================================================"

    #############################################
    # PHASE 1: EXTENDED AMR PROFILING (RGI+CARD)#
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 1: EXTENDED AMR PROFILING (RGI & ABRICATE)      \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m RGI: LOAD CARD DATABASE \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    # Load CARD reference ontology once before per-sample execution
    conda run -n rgi_env rgi load \
    --card_json "${card_db}/card.json" \
    --local

    echo -e "\e[32m CARD database loaded successfully \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    bakta_faa="${ANNOTATION_ECOLI}/bakta/${sample_name}.bakta/${sample_name}.bakta.faa"

    # Checkpoint 1: Validate input assembly
    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] WARNING: Assembly not found: ${filtered_assembly} - skipping \e[0m"
        log_failure "Phase1_AMR" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    RGI_OUT="${RGI_ECOLI}/${sample_name}.rgi"
    ABRICATE_OUT="${ABRICATE_ECOLI}/${sample_name}.abricate"
    mkdir -p "${RGI_OUT}"
    mkdir -p "${ABRICATE_OUT}"

    echo -e "\e[31m ================= \e[0m"
    echo -e "\e[31m RGI: ${sample_name} \e[0m"
    echo -e "\e[31m ================= \e[0m"

    # Contig-based RGI: identifies full resistance genes from assembly contigs
    expected_rgi_contig="${RGI_OUT}/${sample_name}.rgi.txt"
    if [ -s "${expected_rgi_contig}" ]; then
        echo -e "\e[32m   [${sample_name}] Contig RGI already completed — skipping \e[0m"
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
        else
            echo -e "\e[32m Contig RGI complete for ${sample_name} -> ${expected_rgi_contig} \e[0m"
        fi
    fi

    # Protein-based RGI: sensitive detection of point mutations (e.g. gyrA, parC)
    expected_rgi_prot="${RGI_OUT}/${sample_name}.rgi.protein.txt"
    if [ -s "${expected_rgi_prot}" ]; then
        echo -e "\e[32m   [${sample_name}] Protein RGI already completed — skipping \e[0m"
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
        else
            echo -e "\e[32m Protein RGI complete for ${sample_name} -> ${expected_rgi_prot} \e[0m"
        fi
    else
        echo -e "\e[32m   [${sample_name}] Bakta protein FAA not found — skipping protein RGI \e[0m"
        log_failure "Phase1_RGI_protein" "${sample_name}" "INPUT_MISSING" "Bakta protein FAA missing: ${bakta_faa}"
    fi

    echo -e "\e[31m ====================== \e[0m"
    echo -e "\e[31m ABRICATE: ${sample_name} \e[0m"
    echo -e "\e[31m ====================== \e[0m"

    # Rapid screening across curated AMR and virulence databases
    for db in card ncbi resfinder argannot vfdb; do
        abricate_tsv="${ABRICATE_OUT}/${sample_name}.abricate.${db}.tsv"
        if [ -s "${abricate_tsv}" ]; then
            echo -e "\e[32m   [${sample_name}] ABRICATE ${db} already completed — skipping \e[0m"
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
    echo -e "\e[31m ABRICATE: SUMMARIZE ALL E. COLI ISOLATES   \e[0m"
    echo -e "\e[31m =========================================== \e[0m"

    # Aggregate per-database summary tables across all E. coli samples
    for db in card ncbi resfinder argannot vfdb; do
        summary_tsv="${ABRICATE_ECOLI}/all_Ecoli.abricate.${db}.summary.tsv"
        conda run -n BPannotation abricate --summary \
        "${ABRICATE_ECOLI}/"*".abricate/"*".abricate.${db}.tsv" \
        > "${summary_tsv}"
        echo -e "\e[32m ABRICATE summary written: ${summary_tsv} \e[0m"
    done

    echo -e "\e[31m ======================================== \e[0m"
    echo -e "\e[31m RGI: HEATMAP OF ALL E. COLI ISOLATES    \e[0m"
    echo -e "\e[31m ======================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    json_src="${RGI_ECOLI}/${sample_name}.rgi/${sample_name}.rgi.json"
    json_dst="${RGI_ECOLI}/${sample_name}rgi.json"
    if [ -f "${json_src}" ] && [ ! -f "${json_dst}" ]; then
        cp "${json_src}" "${json_dst}"
    fi
done

    # Generate comparative RGI resistance gene presence/absence heatmap
    conda run -n rgi_env rgi heatmap \
    --input "${RGI_ECOLI}" \
    --output "${RGI_ECOLI}/all_Ecoli.rgi_heatmap"
    echo -e "\e[32m RGI heatmap generated: ${RGI_ECOLI}/all_Ecoli.rgi_heatmap \e[0m"

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

    ISESCAN_SAMPLE="${ISESCAN_ECOLI}/${sample_name}.isescan"
    INTEGRON_SAMPLE="${INTEGRON_ECOLI}/${sample_name}.integron"
    ICEFINDER_SAMPLE="${ICEFINDER_ECOLI}/${sample_name}"
    mkdir -p "${ISESCAN_SAMPLE}"
    mkdir -p "${INTEGRON_SAMPLE}"
    mkdir -p "${ICEFINDER_SAMPLE}"

    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m ISESCAN: ${sample_name} \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    # Identify insertion sequences (IS elements, families, terminal repeats)
    if [ -d "${ISESCAN_SAMPLE}/prediction" ]; then
        echo -e "\e[32m   [${sample_name}] ISEScan prediction already exists — skipping \e[0m"
        log_failure "Phase2_ISEScan" "${sample_name}" "SKIPPED_EXISTS" "Directory ${ISESCAN_SAMPLE}/prediction already exists"
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
    if [ -d "${INTEGRON_RESULTS}" ]; then
        echo -e "\e[32m   [${sample_name}] IntegronFinder results already exist — skipping \e[0m"
        log_failure "Phase2_Integron" "${sample_name}" "SKIPPED_EXISTS" "Directory ${INTEGRON_RESULTS} already exists"
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
    cp "${filtered_assembly}" "${ICEFINDER_SAMPLE}/${sample_name}_for_ICEfinder.fasta"
done

    echo -e "\e[32m ICEfinder queries formatted in: ${ICEFINDER_ECOLI}/<sample>/ \e[0m"
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
    renamed_assembly="${PROKKA_ECOLI}/prokka_inputs/${sample_name}.renamed.fasta"
    sample_prokka_dir="${PROKKA_ECOLI}/${sample_name}.prokka"
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
        echo -e "\e[32m   [${sample_name}] Prokka annotation already completed — skipping \e[0m"
        log_failure "Phase3_Prokka" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_prokka_gbk} already exists"
    else
        if ! conda run -n BPannotation prokka \
            --force \
            --cpus "${threads}" \
            --genus Escherichia \
            --species coli \
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

    # Annotate E. coli K-12 MG1655 reference genome if not already annotated
    REF_PROKKA_DIR="${PROKKA_ECOLI}/MG1655.prokka"
    expected_ref_gbk="${REF_PROKKA_DIR}/MG1655.prokka.gbk"
    echo -e "\e[31m =============================== \e[0m"
    echo -e "\e[31m PROKKA: E. COLI K-12 MG1655 REF \e[0m"
    echo -e "\e[31m =============================== \e[0m"

    if [ -s "${Ecoli_ref}" ]; then
        if [ -s "${expected_ref_gbk}" ]; then
            echo -e "\e[32m Reference K-12 MG1655 Prokka GBK already exists — skipping \e[0m"
            log_failure "Phase3_Prokka_Ref" "K12_MG1655" "SKIPPED_EXISTS" "Output ${expected_ref_gbk} already exists"
        else
            if ! conda run -n BPannotation prokka \
                --force \
                --cpus "${threads}" \
                --genus Escherichia \
                --species coli \
                --prefix "MG1655.prokka" \
                --outdir "${REF_PROKKA_DIR}" \
                "${Ecoli_ref}"; then
                echo -e "\e[31m ERROR: Reference Prokka annotation failed \e[0m"
                log_failure "Phase3_Prokka_Ref" "K12_MG1655" "EXECUTION_FAILED" "prokka non-zero exit on reference"
            else
                echo -e "\e[32m Reference K-12 MG1655 Prokka annotation complete -> ${expected_ref_gbk} \e[0m"
            fi
        fi
    else
        echo -e "\e[31m WARNING: E. coli reference genome not found at: ${Ecoli_ref} \e[0m"
        log_failure "Phase3_Prokka_Ref" "K12_MG1655" "INPUT_MISSING" "Reference FASTA not found at ${Ecoli_ref}"
    fi

    #################################################
    # PHASE 4: PATHOGENICITY ISLANDS & PROPHAGES    #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 4: PATHOGENICITY ISLANDS (GIPSY2) & PROPHAGES   \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    # Export LD_LIBRARY_PATH for gipsy2 binary library dependencies
    export LD_LIBRARY_PATH="/storage/student9/miniconda3/envs/gipsy_env/lib:${LD_LIBRARY_PATH:-}"
    Ecoli_ref_gbk="${REF_PROKKA_DIR}/MG1655.prokka.gbk"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    prokka_gbk="${PROKKA_ECOLI}/${sample_name}.prokka/${sample_name}.prokka.gbk"
    gipsy_out="${GIPSY_ECOLI}/${sample_name}.gipsy2"
    phispy_out_dir="${PHISPY_ECOLI}/${sample_name}.phispy"

    if [ ! -s "${prokka_gbk}" ]; then
        echo -e "\e[31m   [${sample_name}] Prokka GBK missing — skipping PAI/prophage analysis \e[0m"
        log_failure "Phase4_PAI" "${sample_name}" "INPUT_MISSING" "Prokka GBK ${prokka_gbk} missing"
        continue
    fi

    mkdir -p "${gipsy_out}"
    mkdir -p "${phispy_out_dir}"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m GIPSY2: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    # Comparative island prediction subtracting E. coli K-12 MG1655 backbone
    expected_gipsy_txt="${gipsy_out}/pathogenicity_islands.txt"
    if [ -s "${Ecoli_ref_gbk}" ]; then
        if [ -s "${expected_gipsy_txt}" ]; then
            echo -e "\e[32m   [${sample_name}] GIPSy2 output already exists — skipping \e[0m"
            log_failure "Phase4_GIPSy2" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_gipsy_txt} already exists"
        else
            if ! conda run -n gipsy_env "${gipsy2}" \
                -q "${prokka_gbk}" \
                -s "${Ecoli_ref_gbk}" \
                -o "${gipsy_out}" \
                -res -vir -met \
                -k fisher \
                --force; then
                echo -e "\e[31m   [${sample_name}] ERROR: GIPSy2 encountered an execution failure \e[0m"
                log_failure "Phase4_GIPSy2" "${sample_name}" "EXECUTION_FAILED" "gipsy2 non-zero exit status"
            else
                echo -e "\e[32m GIPSy2 comparative profiling completed for ${sample_name} \e[0m"
            fi
        fi
    else
        echo -e "\e[31m   [${sample_name}] Reference GBK ${Ecoli_ref_gbk} not available — skipping GIPSy2 \e[0m"
        log_failure "Phase4_GIPSy2" "${sample_name}" "INPUT_MISSING" "Reference GBK ${Ecoli_ref_gbk} missing"
    fi

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PHISPY: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    # Identify prophages integrated in the E. coli genome
    expected_phispy_tsv="${phispy_out_dir}/prophage.tsv"
    if [ -s "${expected_phispy_tsv}" ]; then
        echo -e "\e[32m   [${sample_name}] PhiSpy prediction already exists — skipping \e[0m"
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

    MERGED_PHISPY="${PHISPY_ECOLI}/all_Ecoli.phispy.merged.tsv"
    : > "${MERGED_PHISPY}"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    phispy_tsv="${PHISPY_ECOLI}/${sample_name}.phispy/prophage.tsv"

    if [ -s "${phispy_tsv}" ]; then
        awk -v s="${sample_name}" 'BEGIN{FS=OFS="\t"} {print s, $0}' "${phispy_tsv}" >> "${MERGED_PHISPY}"
        echo -e "\e[32m   ${sample_name}: $(wc -l < "${phispy_tsv}") prophage(s) merged \e[0m"
    fi
done

    echo -e "\e[32m Merged PhiSpy table: ${MERGED_PHISPY} (Total: $(wc -l < "${MERGED_PHISPY}") records) \e[0m"

    #####################################################################
    # PHASE 5: E. COLI SEROTYPE, CLERMONT PHYLOGROUP, MLST & VIRULENCE  #
    #####################################################################

    echo -e "\e[33m ================================================================= \e[0m"
    echo -e "\e[33m PHASE 5: E. COLI TYPING (SEROTYPE, CLERMONT, MLST, VIRULENCE)     \e[0m"
    echo -e "\e[33m ================================================================= \e[0m"

    # 1. Multi-Locus Sequence Typing (MLST) using Achtman 7-gene scheme
    echo -e "\e[31m =================================== \e[0m"
    echo -e "\e[31m MLST (ACHTMAN SCHEME): ALL SAMPLES  \e[0m"
    echo -e "\e[31m =================================== \e[0m"

    ALL_MLST_TSV="${MLST_ECOLI}/all_Ecoli.mlst.tsv"
    if [ -s "${ALL_MLST_TSV}" ]; then
        echo -e "\e[32m MLST summary table already exists — skipping \e[0m"
        log_failure "Phase5_MLST" "ALL_COHORT" "SKIPPED_EXISTS" "Output ${ALL_MLST_TSV} already exists"
    else
        if ! conda run -n BPtyping mlst \
            --scheme ecoli \
            --threads "${threads}" \
            "${CHECKM_INPUTS}/"*.contigs.filtered.fasta \
            > "${ALL_MLST_TSV}"; then
            echo -e "\e[31m ERROR: MLST command failed \e[0m"
            log_failure "Phase5_MLST" "ALL_COHORT" "EXECUTION_FAILED" "mlst command exited non-zero"
        else
            echo -e "\e[32m MLST complete: ${ALL_MLST_TSV} \e[0m"
        fi
    fi

    # 2. O:H Antigen Serotyping (ECTyper)
    # NOTE: Requires environment 'ectyper' with package 'ectyper' installed
    # Install if needed: conda create -n ectyper -c bioconda ectyper
    echo -e "\e[31m =================================== \e[0m"
    echo -e "\e[31m ECTYPER (O:H SEROTYPE): ALL SAMPLES \e[0m"
    echo -e "\e[31m =================================== \e[0m"

    for sample_id in "${sample_ids[@]}"; do
        sample_name="${SAMPLE_PREFIX}${sample_id}"
        filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
        sample_ectyper_out="${ECTYPER_ECOLI}/${sample_name}.ectyper"
        expected_ectyper_tsv="${sample_ectyper_out}/output.tsv"

        if [ ! -s "${filtered_assembly}" ]; then
            echo -e "\e[31m   [${sample_name}] WARNING: input assembly not found — skipping \e[0m"
            log_failure "Phase5_ECTyper" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
            continue
        fi

        if [ -s "${expected_ectyper_tsv}" ]; then
            echo -e "\e[32m   [${sample_name}] ECTyper already completed — skipping \e[0m"
            log_failure "Phase5_ECTyper" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_ectyper_tsv} already exists"
            continue
        fi

        mkdir -p "${sample_ectyper_out}"

        if ! conda run -n ectyper ectyper \
            --input "${filtered_assembly}" \
            --output "${sample_ectyper_out}" \
            --cores "${threads}" \
            --verify; then
            echo -e "\e[31m   [${sample_name}] ERROR: ECTyper execution failed — skipping \e[0m"
            log_failure "Phase5_ECTyper" "${sample_name}" "EXECUTION_FAILED" "ECTyper exited with non-zero status"
            continue
        fi

        if [ ! -s "${expected_ectyper_tsv}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: ECTyper expected output missing or empty \e[0m"
            log_failure "Phase5_ECTyper" "${sample_name}" "OUTPUT_MISSING" "Output ${expected_ectyper_tsv} was not created"
            continue
        fi

        echo -e "\e[32m ECTyper serotyping complete for ${sample_name} -> ${expected_ectyper_tsv} \e[0m"
    done

    # 3. Clermont PCR Phylogrouping (A, B1, B2, C, D, E, F, cryptic clades)
    # NOTE: Requires environment 'clermontyping' with package 'clermontyping' installed
    # Install if needed: conda create -n clermontyping -c bioconda clermontyping
    echo -e "\e[31m ======================================= \e[0m"
    echo -e "\e[31m CLERMONTYPING (PHYLOGROUP): ALL SAMPLES \e[0m"
    echo -e "\e[31m ======================================= \e[0m"

    for sample_id in "${sample_ids[@]}"; do
        sample_name="${SAMPLE_PREFIX}${sample_id}"
        filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
        sample_clermont_out="${CLERMONT_ECOLI}/${sample_name}.clermont"
        expected_clermont_phylogroup="${sample_clermont_out}/phylogroup.txt"

        if [ ! -s "${filtered_assembly}" ]; then
            echo -e "\e[31m   [${sample_name}] WARNING: input assembly not found — skipping \e[0m"
            log_failure "Phase5_Clermontyping" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
            continue
        fi

        if [ -s "${expected_clermont_phylogroup}" ]; then
            echo -e "\e[32m   [${sample_name}] Clermontyping already completed — skipping \e[0m"
            log_failure "Phase5_Clermontyping" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_clermont_phylogroup} already exists"
            continue
        fi

        mkdir -p "${sample_clermont_out}"

        if ! conda run -n clermontyping clermontyping \
            --fasta "${filtered_assembly}" \
            --outdir "${sample_clermont_out}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Clermontyping execution failed — skipping \e[0m"
            log_failure "Phase5_Clermontyping" "${sample_name}" "EXECUTION_FAILED" "Clermontyping exited with non-zero status"
            continue
        fi

        echo -e "\e[32m Clermont phylogrouping complete for ${sample_name} \e[0m"
    done

    # 4. Virulence Determinants (AMRFinderPlus --plus --organism Escherichia)
    echo -e "\e[31m ======================================= \e[0m"
    echo -e "\e[31m AMRFINDERPLUS VIRULENCE: ALL SAMPLES    \e[0m"
    echo -e "\e[31m ======================================= \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    amrfinder_out="${VIRULENCE_ECOLI}/${sample_name}.amrfinder.tsv"

    if [ ! -s "${filtered_assembly}" ]; then
        log_failure "Phase5_AMRFinder" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    if [ -s "${amrfinder_out}" ]; then
        echo -e "\e[32m   [${sample_name}] AMRFinderPlus already completed — skipping \e[0m"
        log_failure "Phase5_AMRFinder" "${sample_name}" "SKIPPED_EXISTS" "Output ${amrfinder_out} already exists"
    else
        if ! conda run -n ncbi amrfinder \
            --nucleotide "${filtered_assembly}" \
            --organism Escherichia \
            --plus \
            --threads "${threads}" \
            --output "${amrfinder_out}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: AMRFinderPlus execution failed \e[0m"
            log_failure "Phase5_AMRFinder" "${sample_name}" "EXECUTION_FAILED" "amrfinder non-zero exit status"
        elif [ ! -s "${amrfinder_out}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: AMRFinderPlus produced empty output \e[0m"
            log_failure "Phase5_AMRFinder" "${sample_name}" "OUTPUT_MISSING" "File ${amrfinder_out} is empty"
        else
            echo -e "\e[32m AMRFinderPlus (virulence+AMR) complete for ${sample_name} -> ${amrfinder_out} \e[0m"
        fi
    fi
done

    echo -e "\e[32m Virulence profiling summary: VFDB summary table generated in Phase 1 \e[0m"
    echo -e "\e[32m   -> ${ABRICATE_ECOLI}/all_Ecoli.abricate.vfdb.summary.tsv \e[0m"

    #################################################
    # PHASE 6: PLASMID IDENTIFICATION & TYPING      #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 6: PLASMID ASSIGNMENT (PLATON, MOB-SUITE, INC)  \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_assembly}" ]; then
        log_failure "Phase6_Plasmid" "${sample_name}" "INPUT_MISSING" "File ${filtered_assembly} missing or empty"
        continue
    fi

    PLATON_OUT="${PLATON_ECOLI}/${sample_name}.platon"
    PLASMIDFINDER_OUT="${PLASMIDFINDER_ECOLI}/${sample_name}.plasmidfinder"
    MOBSUITE_OUT="${MOBSUITE_ECOLI}/${sample_name}.mobsuite"

    mkdir -p "${PLATON_OUT}"
    mkdir -p "${PLASMIDFINDER_OUT}"
    mkdir -p "${MOBSUITE_OUT}"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PLATON: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    # Classify contigs as chromosomal vs plasmid based on marker HMMs and replication proteins
    expected_platon="${PLATON_OUT}/${sample_name}.chromosome.fasta"
    if [ -s "${expected_platon}" ]; then
        echo -e "\e[32m   [${sample_name}] Platon classification already exists — skipping \e[0m"
        log_failure "Phase6_Platon" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_platon} already exists"
    else
        if ! conda run -n plasmid platon \
            --db "${platon_db}" \
            --output "${PLATON_OUT}" \
            --prefix "${sample_name}" \
            --mode sensitivity \
            --threads "${threads}" \
            "${filtered_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Platon execution failed \e[0m"
            log_failure "Phase6_Platon" "${sample_name}" "EXECUTION_FAILED" "platon non-zero exit status"
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
        echo -e "\e[32m   [${sample_name}] PlasmidFinder output already exists — skipping \e[0m"
        log_failure "Phase6_PlasmidFinder" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_plasmidfinder} already exists"
    else
        if ! conda run -n plasmid plasmidfinder.py \
            -i "${filtered_assembly}" \
            -o "${PLASMIDFINDER_OUT}" \
            -p "${plasmidfinder_db}" \
            -l 0.60 \
            -t 0.80 \
            -x; then
            echo -e "\e[31m   [${sample_name}] ERROR: PlasmidFinder execution failed \e[0m"
            log_failure "Phase6_PlasmidFinder" "${sample_name}" "EXECUTION_FAILED" "plasmidfinder.py non-zero exit status"
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
        echo -e "\e[32m   [${sample_name}] MOB-suite results already exist — skipping \e[0m"
        log_failure "Phase6_MOBsuite" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_mob} already exists"
    else
        if ! conda run -n plasmid mob_recon \
            --infile "${filtered_assembly}" \
            --outdir "${MOBSUITE_OUT}" \
            --num_threads "${threads}" \
            --force; then
            echo -e "\e[31m   [${sample_name}] ERROR: MOB-suite execution failed \e[0m"
            log_failure "Phase6_MOBsuite" "${sample_name}" "EXECUTION_FAILED" "mob_recon non-zero exit status"
        else
            echo -e "\e[32m MOB-suite complete for ${sample_name} -> ${expected_mob} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 7: CHIMERA DETECTION & READ COVERAGE    #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 7: READ MAPPING & CHIMERIC BREAKPOINT DETECTION \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    read1t="${PREPROCESSING_ECOLI}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_ECOLI}/${sample_name}.R2.paired.fastq.gz"

    if [ ! -s "${filtered_assembly}" ] || [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        echo -e "\e[31m   [${sample_name}] Fastq or assembly missing — skipping chimera check \e[0m"
        log_failure "Phase7_Chimera" "${sample_name}" "INPUT_MISSING" "Reads or assembly missing for ${sample_name}"
        continue
    fi

    CHIMERA_OUT="${CHIMERA_ECOLI}/${sample_name}.chimera"
    mkdir -p "${CHIMERA_OUT}"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m BWA INDEX: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    conda run -n mapping bwa index "${filtered_assembly}"

    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m BWA MEM: ${sample_name} \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    # Map reads, sort and index BAM
    mapped_bam="${CHIMERA_OUT}/${sample_name}.mapped.bam"
    if [ -s "${mapped_bam}" ]; then
        echo -e "\e[32m   [${sample_name}] BAM already mapped — skipping BWA \e[0m"
        log_failure "Phase7_BWA" "${sample_name}" "SKIPPED_EXISTS" "Output ${mapped_bam} already exists"
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
    conda run -n mapping samtools coverage \
    "${mapped_bam}" \
    > "${CHIMERA_OUT}/${sample_name}.coverage_summary.tsv"

    conda run -n mapping samtools depth \
    -a \
    "${mapped_bam}" \
    > "${CHIMERA_OUT}/${sample_name}.depth_per_base.tsv"

    conda run -n mapping samtools flagstat \
    "${mapped_bam}" \
    > "${CHIMERA_OUT}/${sample_name}.flagstat.txt"

    echo -e "\e[32m Coverage statistics extracted for ${sample_name} \e[0m"

    # Detect coverage breakpoints (misassembly or chimeric joins)
    if [ -f "${detect_coverage_breakpoints}" ]; then
        echo -e "\e[32m Running breakpoint detection for ${sample_name}... \e[0m"
        if ! conda run -n ncbi python3 "${detect_coverage_breakpoints}" \
            --input "${CHIMERA_OUT}/${sample_name}.depth_per_base.tsv" \
            --coverage-summary "${CHIMERA_OUT}/${sample_name}.coverage_summary.tsv" \
            --window 500 \
            --cov-cutoff 0.2 --drop-size 50 \
            --jump-cutoff 3.0 --jump-size 50 \
            --output "${CHIMERA_OUT}/${sample_name}.breakpoints_report.tsv"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Breakpoint detection script failed \e[0m"
            log_failure "Phase7_Breakpoints" "${sample_name}" "EXECUTION_FAILED" "detect_coverage_breakpoints.py error"
        else
            echo -e "\e[32m Breakpoints report generated for ${sample_name} \e[0m"
        fi
    fi
done

    #################################################
    # PHASE 8: RE-ASSEMBLY (UNICYCLER) & QC         #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 8: RE-ASSEMBLY (UNICYCLER) & QUALITY QC (CHECKM)\e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="${SAMPLE_PREFIX}${sample_id}"
    read1t="${PREPROCESSING_ECOLI}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_ECOLI}/${sample_name}.R2.paired.fastq.gz"

    if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        log_failure "Phase8_Unicycler" "${sample_name}" "INPUT_MISSING" "Paired reads missing for ${sample_name}"
        continue
    fi

    UNICYCLER_OUT="${REASSEMBLY_ECOLI}/${sample_name}.unicycler"
    expected_assembly="${UNICYCLER_OUT}/assembly.fasta"
    mkdir -p "${UNICYCLER_OUT}"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m UNICYCLER: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    if [ -s "${expected_assembly}" ]; then
        echo -e "\e[32m   [${sample_name}] Unicycler assembly already exists — skipping \e[0m"
        log_failure "Phase8_Unicycler" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_assembly} already exists"
    else
        if ! conda run -n assembly unicycler \
            -1 "${read1t}" \
            -2 "${read2t}" \
            --mode conservative \
            --threads "${threads}" \
            --out "${UNICYCLER_OUT}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Unicycler execution failed \e[0m"
            log_failure "Phase8_Unicycler" "${sample_name}" "EXECUTION_FAILED" "unicycler exited non-zero"
            continue
        elif [ ! -s "${expected_assembly}" ]; then
            echo -e "\e[31m   [${sample_name}] ERROR: Unicycler expected assembly.fasta missing \e[0m"
            log_failure "Phase8_Unicycler" "${sample_name}" "OUTPUT_MISSING" "File ${expected_assembly} not found"
            continue
        else
            echo -e "\e[32m Unicycler complete for ${sample_name} -> ${expected_assembly} \e[0m"
        fi
    fi

    # QUAST single-isolate assembly metrics
    if [ -s "${expected_assembly}" ]; then
        conda run -n assembly quast.py \
        "${expected_assembly}" \
        --threads "${threads}" \
        --output-dir "${UNICYCLER_OUT}/quast"

        # Copy to multi-sample inputs
        cp -f "${expected_assembly}" \
        "${REASSEMBLY_QUAST_INPUTS_ECOLI}/${sample_name}.unicycler.fasta"
        cp -f "${expected_assembly}" \
        "${REASSEMBLY_ECOLI}/checkm_inputs/${sample_name}.unicycler.fasta"
    fi
done

    # Multi-sample comparative QUAST summary
    echo -e "\e[31m ========================================= \e[0m"
    echo -e "\e[31m QUAST: ALL UNICYCLER ASSEMBLIES (SUMMARY) \e[0m"
    echo -e "\e[31m ========================================= \e[0m"

    if compgen -G "${REASSEMBLY_QUAST_INPUTS_ECOLI}/*.fasta" > /dev/null; then
        conda run -n assembly quast.py \
        --threads "${threads}" \
        --output-dir "${REASSEMBLY_ECOLI}/all_Ecoli.unicycler.quast" \
        "${REASSEMBLY_QUAST_INPUTS_ECOLI}/"*.fasta
        echo -e "\e[32m Multi-sample QUAST complete: ${REASSEMBLY_ECOLI}/all_Ecoli.unicycler.quast \e[0m"
    fi

    # Multi-sample CheckM lineage completeness QC
    echo -e "\e[31m ================================================== \e[0m"
    echo -e "\e[31m CHECKM: ALL UNICYCLER ASSEMBLIES (COMPLETENESS QC)\e[0m"
    echo -e "\e[31m ================================================== \e[0m"

    if compgen -G "${REASSEMBLY_ECOLI}/checkm_inputs/*.fasta" > /dev/null; then
        conda run -n BPstructure checkm lineage_wf \
        -t "${threads}" \
        --reduced_tree \
        --pplacer_threads 1 \
        -x fasta \
        --tab_table \
        "${REASSEMBLY_ECOLI}/checkm_inputs" \
        "${REASSEMBLY_ECOLI}/checkm"

        if [ -f "${REASSEMBLY_ECOLI}/checkm/lineage.ms" ]; then
            conda run -n BPstructure checkm qa \
            "${REASSEMBLY_ECOLI}/checkm/lineage.ms" \
            "${REASSEMBLY_ECOLI}/checkm" \
            -o 2 \
            --tab_table \
            -f "${REASSEMBLY_ECOLI}/checkm/Ecoli.reassembly.quality.checkm.tsv"
            echo -e "\e[32m CheckM QA complete: ${REASSEMBLY_ECOLI}/checkm/Ecoli.reassembly.quality.checkm.tsv \e[0m"
        fi
    fi

    #################################################
    # PIPELINE EXECUTION SUMMARY                    #
    #################################################

    echo ""
    echo -e "\e[32m ================================================================ \e[0m"
    echo -e "\e[32m E. COLI PROFILING & REASSEMBLY COMPLETE — $(date)                \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
    echo ""
    echo -e "\e[32m -- OUTPUT DIRECTORIES ------------------------------------------ \e[0m"
    echo -e "\e[32m  AMR (RGI CARD)         : ${RGI_ECOLI}/ \e[0m"
    echo -e "\e[32m  AMR (ABRICATE multi-db): ${ABRICATE_ECOLI}/ \e[0m"
    echo -e "\e[32m  IS Elements (ISEScan)  : ${ISESCAN_ECOLI}/ \e[0m"
    echo -e "\e[32m  Integrons (IntegronF.) : ${INTEGRON_ECOLI}/ \e[0m"
    echo -e "\e[32m  ICEfinder FASTA inputs : ${ICEFINDER_ECOLI}/ \e[0m"
    echo -e "\e[32m  Prokka Re-annotations  : ${PROKKA_ECOLI}/ \e[0m"
    echo -e "\e[32m  PAIs (GIPSy2)          : ${GIPSY_ECOLI}/ \e[0m"
    echo -e "\e[32m  Prophages (PhiSpy)     : ${PHISPY_ECOLI}/all_Ecoli.phispy.merged.tsv \e[0m"
    echo -e "\e[32m  MLST (Achtman scheme)  : ${MLST_ECOLI}/all_Ecoli.mlst.tsv \e[0m"
    echo -e "\e[32m  Serotype (ECTyper)     : ${ECTYPER_ECOLI}/ \e[0m"
    echo -e "\e[32m  Clermont Phylogroup    : ${CLERMONT_ECOLI}/ \e[0m"
    echo -e "\e[32m  Virulence factors      : ${VIRULENCE_ECOLI}/ \e[0m"
    echo -e "\e[32m  Platon Plasmids        : ${PLATON_ECOLI}/ \e[0m"
    echo -e "\e[32m  PlasmidFinder (Inc)    : ${PLASMIDFINDER_ECOLI}/ \e[0m"
    echo -e "\e[32m  MOB-suite Plasmids     : ${MOBSUITE_ECOLI}/ \e[0m"
    echo -e "\e[32m  Chimera Breakpoint QC  : ${CHIMERA_ECOLI}/ \e[0m"
    echo -e "\e[32m  Unicycler Reassembly   : ${REASSEMBLY_ECOLI}/ \e[0m"
    echo -e "\e[32m  CheckM QC summary      : ${REASSEMBLY_ECOLI}/checkm/Ecoli.reassembly.quality.checkm.tsv \e[0m"
    echo ""
    echo -e "\e[32m -- LOG FILES --------------------------------------------------- \e[0m"
    echo -e "\e[32m  Full execution log     : ${LOG} \e[0m"

    # Report count of failed or skipped samples recorded in FAIL_LOG
    fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING' "${FAIL_LOG}" 2>/dev/null || echo 0)
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
    echo -e "\e[32m  1. If ECTyper or Clermontyping failed due to missing environments: \e[0m"
    echo -e "\e[32m     conda create -n ectyper -c bioconda ectyper \e[0m"
    echo -e "\e[32m     conda create -n clermontyping -c bioconda clermontyping \e[0m"
    echo -e "\e[32m  2. Submit formatted ICEfinder FASTA files to: \e[0m"
    echo -e "\e[32m     https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html \e[0m"
    echo -e "\e[32m  3. Open assembly.gfa in Bandage to confirm circular topology \e[0m"
    echo -e "\e[32m  4. Review breakpoints_report.tsv to inspect any chimeric joins \e[0m"
    echo -e "\e[32m  5. Correlate plasmid replicon types (PlasmidFinder) with MOB-suite \e[0m"
    echo -e "\e[32m     reconstructed plasmids and Platon contig classifications \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"

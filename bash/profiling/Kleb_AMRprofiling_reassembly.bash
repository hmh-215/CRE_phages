#!/bin/bash
# AMR profiling and reassembly (chromosome and plasmids) of K. pneumoniae samples
# Building No. 2
set -euo pipefail

################
# GLOBAL SETUP #
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/Klebsiella_pneumoniae"
PREPROCESSING_KLEB="${WORK_PATH}/preprocessing_Kleb"
CHECKM_INPUTS="${WORK_PATH}/annotation_Kleb/checkm/inputs_Kleb"
ANNOTATION_KLEB="${WORK_PATH}/annotation_Kleb"

threads=16
Kleb_ref="/storage/student9/references/reference_genomes/Klebsiella_pneumoniae/ncbi_dataset/data/GCF_000009885.1/GCF_000009885.1_ASM988v1_genomic.fna"
sample_ids=(01 02 03 04 06 07 10 14 15)

# Paths to tools
gipsy2="/storage/student9/tools/gipsy/gipsy/gipsy2"

# Paths to databases
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

mkdir -p "${RGI_KLEB}"
mkdir -p "${ABRICATE_KLEB}"
mkdir -p "${ISESCAN_KLEB}"
mkdir -p "${INTEGRON_KLEB}"
mkdir -p "${ICEFINDER_KLEB}"
mkdir -p "${PROKKA_KLEB}/prokka_inputs"
mkdir -p "${KLEBORATE_KLEB}"
mkdir -p "${GIPSY_KLEB}"
mkdir -p "${PHISPY_KLEB}"
mkdir -p "${PLATON_KLEB}"
mkdir -p "${PLASMIDFINDER_KLEB}"
mkdir -p "${MOBSUITE_KLEB}"
mkdir -p "${CHIMERA_KLEB}"
mkdir -p "${REASSEMBLY_KLEB}"
mkdir -p "${REASSEMBLY_QUAST_INPUTS_KLEB}"
mkdir -p "${REASSEMBLY_KLEB}/checkm_inputs"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/Kleb_AMRprofiling_reassembly.log"
FAIL_LOG="${WORK_PATH}/Kleb_AMRprofiling_reassembly.failed_skipped.log"
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

echo "================================================================"
echo " K. pneumoniae Extended Profiling & Reassembly Pipeline - Started: $(date)"
echo " Total samples in cohort : ${#sample_ids[@]}"
echo " Reference genome        : ${Kleb_ref}"
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
    sample_name="WS2762512A${sample_id}"
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

    expected_rgi_contig="${RGI_OUT}/${sample_name}.rgi.txt"
    if [ -s "${expected_rgi_contig}" ]; then
        echo -e "\e[32m   [${sample_name}] RGI contig output already exists -> skipping \e[0m"
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
            echo -e "\e[31m   [${sample_name}] ERROR: RGI contig execution failed \e[0m"
            log_failure "Phase1_RGI_contig" "${sample_name}" "EXECUTION_FAILED" "rgi main contig non-zero exit status"
        else
            echo -e "\e[32m RGI contig complete for ${sample_name} \e[0m"
        fi
    fi

    expected_rgi_prot="${RGI_OUT}/${sample_name}.rgi.protein.txt"
    if [ -s "${expected_rgi_prot}" ]; then
        echo -e "\e[32m   [${sample_name}] RGI protein output already exists -> skipping \e[0m"
        log_failure "Phase1_RGI_prot" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_rgi_prot} already exists"
    elif [ ! -s "${bakta_faa}" ]; then
        echo -e "\e[33m   [${sample_name}] Bakta FAA missing (${bakta_faa}) -> skipping protein RGI \e[0m"
        log_failure "Phase1_RGI_prot" "${sample_name}" "INPUT_MISSING" "Bakta FAA missing: ${bakta_faa}"
    else
        if ! conda run -n rgi_env rgi main \
            --input_sequence "${bakta_faa}" \
            --output_file "${RGI_OUT}/${sample_name}.rgi.protein" \
            --input_type protein \
            --alignment_tool BLAST \
            --include_loose \
            --num_threads "${threads}" \
            --clean \
            --local; then
            echo -e "\e[31m   [${sample_name}] ERROR: RGI protein execution failed \e[0m"
            log_failure "Phase1_RGI_prot" "${sample_name}" "EXECUTION_FAILED" "rgi main protein non-zero exit status"
        else
            echo -e "\e[32m RGI protein complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ====================== \e[0m"
    echo -e "\e[31m ABRICATE: ${sample_name} \e[0m"
    echo -e "\e[31m ====================== \e[0m"

    for db in card ncbi resfinder argannot vfdb; do
        expected_abricate="${ABRICATE_OUT}/${sample_name}.abricate.${db}.tsv"
        if [ -s "${expected_abricate}" ]; then
            echo -e "\e[32m   [${sample_name}] ABRICATE ${db} already exists -> skipping \e[0m"
            log_failure "Phase1_ABRICATE_${db}" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_abricate} already exists"
        else
            if ! conda run -n BPannotation abricate \
                --db "${db}" \
                --threads "${threads}" \
                --minid 80 \
                --mincov 80 \
                "${filtered_assembly}" \
                > "${expected_abricate}"; then
                echo -e "\e[31m   [${sample_name}] ERROR: ABRICATE ${db} failed \e[0m"
                log_failure "Phase1_ABRICATE_${db}" "${sample_name}" "EXECUTION_FAILED" "abricate --db ${db} failed"
            else
                echo -e "\e[32m ABRICATE ${db} done for ${sample_name} \e[0m"
            fi
        fi
    done
done

    echo -e "\e[31m ============================================ \e[0m"
    echo -e "\e[31m ABRICATE: SUMMARIZE ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ============================================ \e[0m"

    for db in card ncbi resfinder argannot vfdb; do
        conda run -n BPannotation abricate --summary \
        "${ABRICATE_KLEB}/"*".abricate/"*".abricate.${db}.tsv" \
        > "${ABRICATE_KLEB}/all_Kleb.abricate.${db}.summary.tsv" 2>/dev/null || true
        echo -e "\e[32m ABRICATE summary written: all_Kleb.abricate.${db}.summary.tsv \e[0m"
    done

    echo -e "\e[31m ========================================= \e[0m"
    echo -e "\e[31m RGI: HEATMAP OF ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ========================================= \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="WS2762512A${sample_id}"
    if [ -f "${RGI_KLEB}/${sample_name}.rgi/${sample_name}.rgi.json" ]; then
        mv -f "${RGI_KLEB}/${sample_name}.rgi/${sample_name}.rgi.json" "${RGI_KLEB}/${sample_name}rgi.json"
    fi
done

    if compgen -G "${RGI_KLEB}/*rgi.json" > /dev/null; then
        conda run -n rgi_env rgi heatmap \
        --input "${RGI_KLEB}" \
        --output "${RGI_KLEB}/all_Kleb.rgi_heatmap" 2>/dev/null || true
        echo -e "\e[32m RGI heatmap generated: ${RGI_KLEB}/all_Kleb.rgi_heatmap \e[0m"
    else
        echo -e "\e[33m No RGI JSON files found to generate heatmap \e[0m"
    fi

    #################################################
    # PHASE 2: MOBILE GENETIC ELEMENTS (MGEs)       #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 2: MOBILE GENETIC ELEMENTS (IS, INTEGRONS, ICE) \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="WS2762512A${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Assembly missing -> skipping MGE \e[0m"
        log_failure "Phase2_MGE" "${sample_name}" "INPUT_MISSING" "Assembly ${filtered_assembly} missing"
        continue
    fi

    mkdir -p "${ISESCAN_KLEB}/${sample_name}.isescan"
    mkdir -p "${INTEGRON_KLEB}/${sample_name}.integron"
    mkdir -p "${ICEFINDER_KLEB}/${sample_name}"

    ISESCAN_KLEB_SAMPLE="${ISESCAN_KLEB}/${sample_name}.isescan"
    INTEGRON_KLEB_SAMPLE="${INTEGRON_KLEB}/${sample_name}.integron"

    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m ISESCAN: ${sample_name} \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    expected_isescan="${ISESCAN_KLEB_SAMPLE}/${sample_name}.contigs.filtered.fasta.csv"
    if [ -s "${expected_isescan}" ]; then
        echo -e "\e[32m   [${sample_name}] ISEScan output already exists -> skipping \e[0m"
        log_failure "Phase2_ISEScan" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_isescan} already exists"
    else
        if ! conda run -n recombination isescan.py \
            --nthread "${threads}" \
            --seqfile "${filtered_assembly}" \
            --output "${ISESCAN_KLEB_SAMPLE}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: ISEScan execution failed \e[0m"
            log_failure "Phase2_ISEScan" "${sample_name}" "EXECUTION_FAILED" "isescan.py non-zero exit status"
        else
            echo -e "\e[32m ISEScan complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ============================= \e[0m"
    echo -e "\e[31m INTEGRON FINDER: ${sample_name} \e[0m"
    echo -e "\e[31m ============================= \e[0m"

    INTEGRON_RESULTS="${INTEGRON_KLEB_SAMPLE}/Results_Integron_Finder_${sample_name}.contigs.filtered"
    if [ -d "${INTEGRON_RESULTS}" ] && [ "$(ls -A "${INTEGRON_RESULTS}" 2>/dev/null)" ]; then
        echo -e "\e[32m   [${sample_name}] IntegronFinder results already exist -> skipping \e[0m"
        log_failure "Phase2_IntegronFinder" "${sample_name}" "SKIPPED_EXISTS" "Directory ${INTEGRON_RESULTS} already populated"
    else
        if ! conda run -n recombination integron_finder \
            --gbk --pdf \
            --circ \
            --outdir "${INTEGRON_KLEB_SAMPLE}" \
            "${filtered_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: IntegronFinder execution failed \e[0m"
            log_failure "Phase2_IntegronFinder" "${sample_name}" "EXECUTION_FAILED" "integron_finder non-zero exit status"
        else
            if [ -d "${INTEGRON_RESULTS}" ]; then
                for contig_gbk in "${INTEGRON_RESULTS}"/contig_*.gbk; do
                    [ -f "${contig_gbk}" ] || continue
                    contig_name=$(basename "${contig_gbk}" .gbk)
                    mv "${contig_gbk}" "${INTEGRON_RESULTS}/${sample_name}.${contig_name}.integron.gbk"
                done
                pdf_count=0
                for pdf_file in "${INTEGRON_RESULTS}"/contig_*_*.pdf; do
                    [ -f "${pdf_file}" ] || continue
                    integron_num=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f1 | rev)
                    contig_id=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f2- | rev | sed 's/^contig_//')
                    new_name="${sample_name}.contig${contig_id}.integron_${integron_num}.pdf"
                    mv "${pdf_file}" "${INTEGRON_RESULTS}/${new_name}"
                    pdf_count=$(( pdf_count + 1 ))
                done
                echo -e "\e[32m IntegronFinder PDFs renamed: ${pdf_count} \e[0m"
            fi
            echo -e "\e[32m IntegronFinder complete for ${sample_name} \e[0m"
        fi
    fi

    # Prepare ICEfinder submission file
    expected_ice="${ICEFINDER_KLEB}/${sample_name}/${sample_name}_for_ICEfinder.fasta"
    if [ -s "${expected_ice}" ]; then
        echo -e "\e[32m   [${sample_name}] ICEfinder input already formatted -> skipping \e[0m"
        log_failure "Phase2_ICEfinder" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_ice} already exists"
    else
        cp "${filtered_assembly}" "${expected_ice}"
        echo -e "\e[32m ICEfinder submission file ready: ${expected_ice} \e[0m"
    fi
done

    #############################################
    # PHASE 3: GENOME RE-ANNOTATION (PROKKA)    #
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 3: GENOME RE-ANNOTATION (PROKKA)                 \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="WS2762512A${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    renamed_assembly="${PROKKA_KLEB}/prokka_inputs/${sample_name}.renamed.fasta"
    sample_prokka_dir="${PROKKA_KLEB}/${sample_name}.prokka"
    expected_prokka_gff="${sample_prokka_dir}/${sample_name}.prokka.gff"

    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Assembly missing -> skipping Prokka \e[0m"
        log_failure "Phase3_Prokka" "${sample_name}" "INPUT_MISSING" "Assembly ${filtered_assembly} missing"
        continue
    fi

    if [ -s "${expected_prokka_gff}" ]; then
        echo -e "\e[32m   [${sample_name}] Prokka output already exists -> skipping \e[0m"
        log_failure "Phase3_Prokka" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_prokka_gff} already exists"
    else
        awk '/^>/ {
            counter++
            print ">contig_" counter
            next
        }
        { print }' "${filtered_assembly}" > "${renamed_assembly}"

        echo -e "\e[31m ==================== \e[0m"
        echo -e "\e[31m PROKKA: ${sample_name} \e[0m"
        echo -e "\e[31m ==================== \e[0m"

        if ! conda run -n BPannotation prokka \
            --force \
            --cpus "${threads}" \
            --genus Klebsiella \
            --species pneumoniae \
            --prefix "${sample_name}.prokka" \
            --outdir "${sample_prokka_dir}" \
            "${renamed_assembly}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Prokka failed \e[0m"
            log_failure "Phase3_Prokka" "${sample_name}" "EXECUTION_FAILED" "prokka non-zero exit status"
        else
            echo -e "\e[32m Prokka complete for ${sample_name} \e[0m"
        fi
    fi
done

    echo -e "\e[31m ============================ \e[0m"
    echo -e "\e[31m PROKKA: REFERENCE NTUH-K2004 \e[0m"
    echo -e "\e[31m ============================ \e[0m"

    expected_ref_gff="${PROKKA_KLEB}/NTUH-K2004.prokka/NTUH-K2004.prokka.gff"
    if [ -s "${expected_ref_gff}" ]; then
        echo -e "\e[32m Reference NTUH-K2004 Prokka annotation already exists -> skipping \e[0m"
        log_failure "Phase3_Prokka_Ref" "NTUH-K2004" "SKIPPED_EXISTS" "Output ${expected_ref_gff} already exists"
    elif [ ! -s "${Kleb_ref}" ]; then
        echo -e "\e[31m Reference FASTA missing: ${Kleb_ref} -> skipping reference Prokka \e[0m"
        log_failure "Phase3_Prokka_Ref" "NTUH-K2004" "INPUT_MISSING" "Reference FASTA missing: ${Kleb_ref}"
    else
        conda run -n BPannotation prokka \
        --force \
        --cpus "${threads}" \
        --genus Klebsiella \
        --species pneumoniae \
        --prefix NTUH-K2004.prokka \
        --outdir "${PROKKA_KLEB}/NTUH-K2004.prokka" \
        "${Kleb_ref}"
        echo -e "\e[32m Prokka complete for NTUH-K2004 \e[0m"
    fi

    #################################################
    # PHASE 4: PATHOGENICITY ISLANDS & VIRULENCE    #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 4: PATHOGENICITY ISLANDS (GIPSY2) & PROPHAGES   \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="WS2762512A${sample_id}"
    prokka_gbk="${PROKKA_KLEB}/${sample_name}.prokka/${sample_name}.prokka.gbk"
    Kleb_ref_gbk="${PROKKA_KLEB}/NTUH-K2004.prokka/NTUH-K2004.prokka.gbk"

    mkdir -p "${GIPSY_KLEB}/${sample_name}.gipsy2"
    mkdir -p "${PHISPY_KLEB}/${sample_name}.phispy"
    PHISPY_OUT="${PHISPY_KLEB}/${sample_name}.phispy"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m GIPSY2: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    expected_gipsy="${GIPSY_KLEB}/${sample_name}.gipsy2"
    if [ ! -s "${prokka_gbk}" ] || [ ! -s "${Kleb_ref_gbk}" ]; then
        echo -e "\e[31m   [${sample_name}] Prokka GBK missing for sample or reference -> skipping GIPSy2 \e[0m"
        log_failure "Phase4_GIPSy2" "${sample_name}" "INPUT_MISSING" "Prokka GBK missing for ${sample_name} or reference"
    elif [ -d "${expected_gipsy}" ] && [ "$(ls -A "${expected_gipsy}" 2>/dev/null)" ]; then
        echo -e "\e[32m   [${sample_name}] GIPSy2 results already exist -> skipping \e[0m"
        log_failure "Phase4_GIPSy2" "${sample_name}" "SKIPPED_EXISTS" "Directory ${expected_gipsy} already populated"
    else
        export LD_LIBRARY_PATH="/storage/student9/miniconda3/envs/gipsy_env/lib:${LD_LIBRARY_PATH:-}"
        if ! conda run -n gipsy_env "${gipsy2}" \
            -q "${prokka_gbk}" \
            -s "${Kleb_ref_gbk}" \
            -o "${GIPSY_KLEB}/${sample_name}.gipsy2" \
            -res -vir -met \
            -k fisher \
            --force; then
            echo -e "\e[31m   [${sample_name}] ERROR: GIPSy2 execution failed \e[0m"
            log_failure "Phase4_GIPSy2" "${sample_name}" "EXECUTION_FAILED" "gipsy2 non-zero exit status"
        else
            echo -e "\e[32m GIPSY2 complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PHISPY: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    expected_phispy="${PHISPY_OUT}/prophage_coordinates.tsv"
    if [ -s "${expected_phispy}" ]; then
        echo -e "\e[32m   [${sample_name}] PhiSpy output already exists -> skipping \e[0m"
        log_failure "Phase4_PhiSpy" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_phispy} already exists"
    elif [ ! -s "${prokka_gbk}" ]; then
        echo -e "\e[31m   [${sample_name}] Prokka GBK missing -> skipping PhiSpy \e[0m"
        log_failure "Phase4_PhiSpy" "${sample_name}" "INPUT_MISSING" "Prokka GBK missing: ${prokka_gbk}"
    else
        if ! conda run -n recombination PhiSpy.py \
            "${prokka_gbk}" \
            -o "${PHISPY_OUT}" \
            --output_choice 7 \
            --threads "${threads}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: PhiSpy execution failed \e[0m"
            log_failure "Phase4_PhiSpy" "${sample_name}" "EXECUTION_FAILED" "PhiSpy.py non-zero exit status"
        else
            echo -e "\e[32m PhiSpy complete for ${sample_name} \e[0m"
        fi
    fi
done

    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m KLEBORATE: ALL K. PNEUMONIAE SAMPLES \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    expected_kleborate="${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv"
    if [ -s "${expected_kleborate}" ]; then
        echo -e "\e[32m Kleborate output already exists -> skipping \e[0m"
        log_failure "Phase4_Kleborate" "ALL" "SKIPPED_EXISTS" "Output ${expected_kleborate} already exists"
    else
        conda run -n BPtyping kleborate \
        --threads "${threads}" \
        --preset kpsc \
        --assemblies "${CHECKM_INPUTS}/"*.fasta \
        --outfile "${expected_kleborate}" || true
        echo -e "\e[32m Kleborate complete: ${expected_kleborate} \e[0m"
    fi

    #################################################
    # PHASE 5: PLASMID IDENTIFICATION & TYPING      #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 5: PLASMID ASSIGNMENT (PLATON, MOB-SUITE, INC)  \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
    sample_name="WS2762512A${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"

    if [ ! -s "${filtered_assembly}" ]; then
        echo -e "\e[31m   [${sample_name}] Assembly missing -> skipping plasmid assignment \e[0m"
        log_failure "Phase5_Plasmids" "${sample_name}" "INPUT_MISSING" "Assembly ${filtered_assembly} missing"
        continue
    fi

    PLATON_OUT="${PLATON_KLEB}/${sample_name}.platon"
    PLASMIDFINDER_OUT="${PLASMIDFINDER_KLEB}/${sample_name}.plasmidfinder"
    MOBSUITE_OUT="${MOBSUITE_KLEB}/${sample_name}.mobsuite"

    mkdir -p "${PLATON_OUT}"
    mkdir -p "${PLASMIDFINDER_OUT}"
    mkdir -p "${MOBSUITE_OUT}"

    echo -e "\e[31m ==================== \e[0m"
    echo -e "\e[31m PLATON: ${sample_name} \e[0m"
    echo -e "\e[31m ==================== \e[0m"

    expected_platon="${PLATON_OUT}/${sample_name}.tsv"
    if [ -s "${expected_platon}" ]; then
        echo -e "\e[32m   [${sample_name}] Platon output already exists -> skipping \e[0m"
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

    expected_plasmidfinder="${PLASMIDFINDER_OUT}/results_tab.tsv"
    if [ -s "${expected_plasmidfinder}" ]; then
        echo -e "\e[32m   [${sample_name}] PlasmidFinder output already exists -> skipping \e[0m"
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
            echo -e "\e[32m PlasmidFinder complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m MOB-SUITE: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    expected_mob="${MOBSUITE_OUT}/mobtyper_results.txt"
    if [ -s "${expected_mob}" ]; then
        echo -e "\e[32m   [${sample_name}] MOB-suite results already exist -> skipping \e[0m"
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
            echo -e "\e[32m MOB-suite complete for ${sample_name} \e[0m"
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
    sample_name="WS2762512A${sample_id}"
    filtered_assembly="${CHECKM_INPUTS}/${sample_name}.contigs.filtered.fasta"
    read1t="${PREPROCESSING_KLEB}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_KLEB}/${sample_name}.R2.paired.fastq.gz"

    if [ ! -s "${filtered_assembly}" ] || [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        echo -e "\e[31m   [${sample_name}] Fastq or assembly missing -> skipping chimera check \e[0m"
        log_failure "Phase6_Chimera" "${sample_name}" "INPUT_MISSING" "Reads or assembly missing for ${sample_name}"
        continue
    fi

    CHIMERA_OUT="${CHIMERA_KLEB}/${sample_name}.chimera"
    mkdir -p "${CHIMERA_OUT}"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m BWA INDEX: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    conda run -n mapping bwa index "${filtered_assembly}"

    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m BWA MEM: ${sample_name} \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    mapped_bam="${CHIMERA_OUT}/${sample_name}.mapped.bam"
    if [ -s "${mapped_bam}" ]; then
        echo -e "\e[32m   [${sample_name}] BAM already mapped -> skipping BWA \e[0m"
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

    echo -e "\e[32m Coverage stats written for ${sample_name} \e[0m"
done

    echo -e "\e[32m Mapping rate summary (all samples): \e[0m"
    for sample_id in "${sample_ids[@]}"; do
        sample_name="WS2762512A${sample_id}"
        flagstat_file="${CHIMERA_KLEB}/${sample_name}.chimera/${sample_name}.flagstat.txt"
        if [ -s "${flagstat_file}" ]; then
            rate=$(grep "mapped (" "${flagstat_file}" | head -1 | awk '{print $5}' | tr -d '()')
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
    sample_name="WS2762512A${sample_id}"
    read1t="${PREPROCESSING_KLEB}/${sample_name}.R1.paired.fastq.gz"
    read2t="${PREPROCESSING_KLEB}/${sample_name}.R2.paired.fastq.gz"

    if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
        echo -e "\e[31m   [${sample_name}] Trimmed reads missing -> skipping Unicycler \e[0m"
        log_failure "Phase7_Unicycler" "${sample_name}" "INPUT_MISSING" "Reads missing for ${sample_name}"
        continue
    fi

    UNICYCLER_OUT="${REASSEMBLY_KLEB}/${sample_name}.unicycler"
    mkdir -p "${UNICYCLER_OUT}"

    echo -e "\e[31m ======================= \e[0m"
    echo -e "\e[31m UNICYCLER: ${sample_name} \e[0m"
    echo -e "\e[31m ======================= \e[0m"

    expected_unicycler="${UNICYCLER_OUT}/assembly.fasta"
    if [ -s "${expected_unicycler}" ]; then
        echo -e "\e[32m   [${sample_name}] Unicycler assembly already exists -> skipping \e[0m"
        log_failure "Phase7_Unicycler" "${sample_name}" "SKIPPED_EXISTS" "Output ${expected_unicycler} already exists"
    else
        if ! conda run -n assembly unicycler \
            -1 "${read1t}" \
            -2 "${read2t}" \
            --mode conservative \
            --threads "${threads}" \
            --out "${UNICYCLER_OUT}"; then
            echo -e "\e[31m   [${sample_name}] ERROR: Unicycler execution failed \e[0m"
            log_failure "Phase7_Unicycler" "${sample_name}" "EXECUTION_FAILED" "unicycler non-zero exit status"
        else
            echo -e "\e[32m Unicycler complete for ${sample_name} \e[0m"
        fi
    fi

    echo -e "\e[31m =================== \e[0m"
    echo -e "\e[31m QUAST: ${sample_name} \e[0m"
    echo -e "\e[31m =================== \e[0m"

    if [ -s "${expected_unicycler}" ]; then
        conda run -n assembly quast.py \
        "${expected_unicycler}" \
        --threads "${threads}" \
        --output-dir "${UNICYCLER_OUT}/quast" 2>/dev/null || true
        echo -e "\e[32m QUAST complete for ${sample_name} \e[0m"
    fi
done

for sample_id in "${sample_ids[@]}"; do
    sample_name="WS2762512A${sample_id}"
    if [ -s "${REASSEMBLY_KLEB}/${sample_name}.unicycler/assembly.fasta" ]; then
        cp -f "${REASSEMBLY_KLEB}/${sample_name}.unicycler/assembly.fasta" \
        "${REASSEMBLY_QUAST_INPUTS_KLEB}/${sample_name}.unicycler.fasta"
        cp -f "${REASSEMBLY_KLEB}/${sample_name}.unicycler/assembly.fasta" \
        "${REASSEMBLY_KLEB}/checkm_inputs/${sample_name}.unicycler.fasta"
    fi
done

    echo -e "\e[31m ========================================= \e[0m"
    echo -e "\e[31m QUAST: ALL UNICYCLER ASSEMBLIES (SUMMARY) \e[0m"
    echo -e "\e[31m ========================================= \e[0m"

    if compgen -G "${REASSEMBLY_QUAST_INPUTS_KLEB}/*.fasta" > /dev/null; then
        conda run -n assembly quast.py \
        --threads "${threads}" \
        --output-dir "${REASSEMBLY_KLEB}/all_Kleb.unicycler.quast" \
        "${REASSEMBLY_QUAST_INPUTS_KLEB}/"*.fasta 2>/dev/null || true
        echo -e "\e[32m Multi-sample QUAST summary written \e[0m"
    fi

    echo -e "\e[31m ================================================== \e[0m"
    echo -e "\e[31m CHECKM: ALL UNICYCLER ASSEMBLIES (COMPLETENESS QC) \e[0m"
    echo -e "\e[31m ================================================== \e[0m"

    if compgen -G "${REASSEMBLY_KLEB}/checkm_inputs/*.fasta" > /dev/null; then
        conda run -n BPstructure checkm lineage_wf \
        -t "${threads}" \
        --reduced_tree \
        --pplacer_threads 1 \
        -x fasta \
        --tab_table \
        "${REASSEMBLY_KLEB}/checkm_inputs" \
        "${REASSEMBLY_KLEB}/checkm" 2>/dev/null || true

        if [ -f "${REASSEMBLY_KLEB}/checkm/lineage.ms" ]; then
            conda run -n BPstructure checkm qa \
            "${REASSEMBLY_KLEB}/checkm/lineage.ms" \
            "${REASSEMBLY_KLEB}/checkm" \
            -o 2 \
            --tab_table \
            -f "${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv"
            echo -e "\e[32m CheckM complete: ${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv \e[0m"
        else
            echo -e "\e[33m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
        fi
    fi

    #################################################
    # PIPELINE EXECUTION SUMMARY                    #
    #################################################

    echo ""
    echo -e "\e[32m ================================================================ \e[0m"
    echo -e "\e[32m K. PNEUMONIAE PROFILING & REASSEMBLY COMPLETE â€” $(date)          \e[0m"
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
    echo -e "\e[32m  Prophages (PhiSpy)     : ${PHISPY_KLEB}/ \e[0m"
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
    echo -e "\e[32m  1. Submit formatted ICEfinder FASTA files to: \e[0m"
    echo -e "\e[32m     https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html \e[0m"
    echo -e "\e[32m  2. Open assembly.gfa in Bandage to confirm circular topology \e[0m"
    echo -e "\e[32m  3. Review depth_per_base.tsv to inspect any chimeric joins \e[0m"
    echo -e "\e[32m  4. Correlate plasmid replicon types (PlasmidFinder) with MOB-suite \e[0m"
    echo -e "\e[32m     reconstructed plasmids and Platon contig classifications \e[0m"
    echo -e "\e[32m  5. Check Kleborate output for yersiniabactin/colibactin/aerobactin \e[0m"
    echo -e "\e[32m ================================================================ \e[0m"
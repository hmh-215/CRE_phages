#!/bin/bash
# Merge all Thai Nguyen (TN) hospital E. coli datasets across batch 2 and batch 3
# Workflow (collate_*.py edition):
#   Phase 0 : Cohort list + per-sample manifests (sample<TAB>path) + stage Prokka .gff files
#             (Prokka only: Panaroo crashes on Bakta .gff3)
#   Phase 1 : AMR collation        -> collate_amr.py (ABRicate CARD + AMRFinderPlus)
#   Phase 2 : Panaroo pangenome    -> panaroo_Ecoli_TN/results/
#   Phase 3 : Annotation collation -> collate_annotations.py
#             - Step 3.1: Bakta GFF3 + Bakta eggNOG        (GFF mode)
#             - Step 3.2: Panaroo + Prokka GFF + Prokka eggNOG (Panaroo mode)
#   Phase 4 : Genome overview      -> collate_genome_overview.py
#             (QUAST + CheckM + Panaroo + MLST + ECTyper + ISEScan/IntegronFinder/PhiSpy)
# Replaces: merge_abricate / merge_amrfinder / merge_abricate_amrfinder / merge_eggnog /
#           merge_quast / merge_mge / merge_genome_overview / merge_pangenome_annotations
#           and the CheckM head/tail concatenation.
# Building no. 9
set -euo pipefail

################
# GLOBAL SETUP #
################

# Base project locations (supports either bacteria_phages_X or bacteria_phage_X naming)
BP2_PATH="/storage/student9/projects/bacteria_phages_2"
[ ! -d "${BP2_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage_2" ] && \
    BP2_PATH="/storage/student9/projects/bacteria_phage_2"

BP3_PATH="/storage/student9/projects/bacteria_phages_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage_3" ] && \
    BP3_PATH="/storage/student9/projects/bacteria_phage_3"

TOOL_PATH="/storage/student9/tools"

# Target reorganization root folder
MERGE_PATH="${BP3_PATH}/merge_Ecoli_TN"

# Subfolder organization:
# merge_Ecoli_TN/
#   |___ paths/                    (cohort list + manifests)
#   |___ panaroo_Ecoli_TN/
#   |     |___ gff/
#   |     |___ results/
#   |___ merge_individual_outputs/ (per-tool matrices, eggNOG summaries, overview details)
#   |___ merge_metadata/           (final collated tables)
PATHS_DIR="${MERGE_PATH}/paths"
PANAROO_DIR="${MERGE_PATH}/panaroo_Ecoli_TN"
PANAROO_GFF_DIR="${PANAROO_DIR}/gff"
PANAROO_RES_DIR="${PANAROO_DIR}/results"
INDIVIDUAL_OUT_DIR="${MERGE_PATH}/merge_individual_outputs"
METADATA_OUT_DIR="${MERGE_PATH}/merge_metadata"

mkdir -p "${MERGE_PATH}"
mkdir -p "${PATHS_DIR}"
mkdir -p "${PANAROO_GFF_DIR}"
mkdir -p "${PANAROO_RES_DIR}"
mkdir -p "${INDIVIDUAL_OUT_DIR}"
mkdir -p "${METADATA_OUT_DIR}"

# Cohort: isolates of each batch (order is preserved in the collated tables)
sample_ids_2=(
    132895-LJF30331 134194-LJF30332 144342-LJF30333 156589-LJF30334 194083-LJF30335
    216381-LJF30336 216607-LJF30337 243686-LJF30338 252837-LJF30339 257625-LJF30340
    319707-LJF30341 332413-LJF30342 346152-LJF30343 125919-LJF30344 131695-LJF30345
    196951-LJF30346 222949-LJF30347 241208-LJF30348 251713-LJF30349 247233-LJF30350
    256139-LJF30351 258569-LJF30352
)
sample_ids_3=(
    WS2762607A01 WS2762607A02 WS2762607A03 WS2762607A04
    WS2762607A06 WS2762607A07 WS2762607A08 WS2762607A09 WS2762607A10
    WS2762607A11 WS2762607A12 WS2762607A13 WS2762607A14 WS2762607A15
    WS2762607A16 WS2762607A17 WS2762607A18 WS2762607A19 WS2762607A20
    WS2762607A21 WS2762607A22 WS2762607A23 WS2762607A24 WS2762607A25
    WS2762607A26 WS2762607A27 WS2762607A28 WS2762607A29 WS2762607A30
    WS2762607A31 WS2762607A32 WS2762607A33 WS2762607A34 WS2762607A35
    WS2762607A36 WS2762607A37 WS2762607A38 WS2762607A39 WS2762607A40
    WS2762607A41 WS2762607A42 WS2762607A43 WS2762607A44 WS2762607A45
    WS2762607A46 WS2762607A47 WS2762607A48 WS2762607A53
)

# Cohort list + manifests (Phase 0 writes them; every collate_*.py reads them)
# Manifest format: one "sample<TAB>path" per line; the cohort file is one sample per line.
samples_list="${PATHS_DIR}/samples_TN.txt"
paths_abricate="${PATHS_DIR}/paths_abricate_card_TN.txt"
paths_amrfinder="${PATHS_DIR}/paths_amrfinder_TN.txt"
paths_eggnog_bakta="${PATHS_DIR}/paths_eggnog_bakta_TN.txt"
paths_eggnog_prokka="${PATHS_DIR}/paths_eggnog_prokka_TN.txt"
paths_bakta_gff="${PATHS_DIR}/paths_bakta_gff_TN.txt"
paths_quast="${PATHS_DIR}/paths_quast_TN.txt"
paths_isescan="${PATHS_DIR}/paths_isescan_TN.txt"
paths_integron="${PATHS_DIR}/paths_integron_TN.txt"
paths_phispy="${PATHS_DIR}/paths_phispy_TN.txt"
paths_ectyper="${PATHS_DIR}/paths_ectyper_TN.txt"
paths_platon="${PATHS_DIR}/paths_platon_TN.txt"
paths_plasmidfinder="${PATHS_DIR}/paths_plasmidfinder_TN.txt"
paths_mobsuite="${PATHS_DIR}/paths_mobsuite_TN.txt"
paths_breakpoints="${PATHS_DIR}/paths_breakpoints_TN.txt"
paths_annotation_summary="${PATHS_DIR}/paths_annotation_summary_TN.txt"

# Batch-wide tables passed straight to collate_genome_overview.py (no manifest needed)
checkm_bp2="${BP2_PATH}/Escherichia_coli/annotation_Ecoli/checkm/all_Ecoli.quality.checkm.tsv"
checkm_bp3="${BP3_PATH}/Escherichia_coli/annotation_Ecoli/checkm/all_Ecoli.quality.checkm.tsv"
mlst_bp2="${BP2_PATH}/Escherichia_coli/typing_Ecoli/mlst/all_Ecoli.mlst.tsv"
mlst_bp3="${BP3_PATH}/Escherichia_coli/typing_Ecoli/mlst/all_Ecoli.mlst.tsv"

# Output files
panaroo_matrix="${PANAROO_RES_DIR}/gene_presence_absence.csv"
out_amr="${METADATA_OUT_DIR}/amr_reconciled_TN.tsv"
amr_per_tool_prefix="${INDIVIDUAL_OUT_DIR}/amr_TN"
out_bakta_eggnog="${METADATA_OUT_DIR}/bakta_eggnog_annotated_genes_TN.tsv"
eggnog_bakta_prefix="${INDIVIDUAL_OUT_DIR}/eggnog_bakta"
out_panaroo_prokka="${METADATA_OUT_DIR}/panaroo_prokka_annotated_genes_TN.tsv"
out_cluster_summary="${METADATA_OUT_DIR}/panaroo_prokka_clusters_summary_TN.tsv"
eggnog_prokka_prefix="${INDIVIDUAL_OUT_DIR}/eggnog_prokka"
out_overview="${METADATA_OUT_DIR}/genome_overview_qc_TN.tsv"
overview_detail_prefix="${INDIVIDUAL_OUT_DIR}/overview_TN"

# Tunable parameters
threads=16                  # Panaroo threads
force_collate=false         # true = regenerate collated tables even if they already exist
abricate_database="card"    # ABRicate database label kept from the ABRicate manifest
amrfinder_type_filter=()    # e.g. (AMR) keeps only Type==AMR rows; empty = keep all types
core_cutoff=0.99            # Panaroo compartments: Core    >= this cohort frequency
soft_core_cutoff=0.95       #                       SoftCore >= this
shell_cutoff=0.15           #                       Shell    >= this (below = Cloud)
species="ecoli"             # overview size range: ecoli | klebsiella | pseudomonas
min_completeness=95.0       # CheckM PASS threshold (%)
max_contamination=5.0       # CheckM PASS threshold (%)
max_contigs=300             # Assembly_Status PASS threshold

# Custom Python helper scripts (all run in the 'ncbi' env; requires Python >= 3.10)
collate_amr="${TOOL_PATH}/collate_amr.py"
collate_annotations="${TOOL_PATH}/collate_annotations.py"
collate_genome_overview="${TOOL_PATH}/collate_genome_overview.py"

# Dual-log setup
LOG="${MERGE_PATH}/merge_script_TN.log"
FAIL_LOG="${MERGE_PATH}/merge_script_TN.failed_skipped.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================" >> "${FAIL_LOG}"
echo " Failure & Skip Log — Started: $(date)" >> "${FAIL_LOG}"
echo "======================================================" >> "${FAIL_LOG}"

# Helper function to record failures/skips to FAIL_LOG with timestamp
log_failure() {
    local phase="$1"
    local sample="$2"
    local status="$3"  # e.g. "SKIPPED_EXISTS", "INPUT_MISSING", "EXECUTION_FAILED", "OUTPUT_MISSING"
    local reason="$4"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${phase}] [${sample}] [${status}] ${reason}" >> "${FAIL_LOG}"
}

echo "======================================================"
echo " Merge all Thai Nguyen hospital E. coli samples — Started: $(date)"
echo " Cohort         : ${#sample_ids_2[@]} (batch 2) + ${#sample_ids_3[@]} (batch 3) isolates"
echo " Root directory : ${MERGE_PATH}"
echo " Collate scripts: ${TOOL_PATH}/collate_{amr,annotations,genome_overview}.py"
echo " Full log       : ${LOG}"
echo " Failure log    : ${FAIL_LOG}"
echo "======================================================"

########################################################
# PHASE 0: COHORT LIST, MANIFESTS & PROKKA GFF STAGING #
########################################################

echo -e "\e[33m ====================================================== \e[0m"
echo -e "\e[33m PHASE 0: COHORT LIST, MANIFESTS & GFF STAGING            \e[0m"
echo -e "\e[33m ====================================================== \e[0m"

# add_entry <manifest> <phase_tag> <sample> <path>
# Appends "sample<TAB>path". The entry is kept even if the file is missing so that the
# collate scripts report "not run" (NA / empty cell) instead of "none found"; the missing
# input is also written to FAIL_LOG.
add_entry() {
    local manifest="$1" tag="$2" sample="$3" path="$4"
    printf '%s\t%s\n' "${sample}" "${path}" >> "${manifest}"
    if [ ! -e "${path}" ]; then
        log_failure "${tag}" "${sample}" "INPUT_MISSING" "Expected input not found: ${path}"
    fi
}

# manifest_has_existing <manifest>: succeeds if at least one listed path exists
manifest_has_existing() {
    local manifest="$1" _sample _path
    [ -s "${manifest}" ] || return 1
    while IFS=$'\t' read -r _sample _path; do
        if [ -n "${_path}" ] && [ -e "${_path}" ]; then
            return 0
        fi
    done < "${manifest}"
    return 1
}

# stage_batch <label> <batch_root> <sample_id>...
# Records every per-sample path in the manifests and copies the Prokka GFF for Panaroo.
stage_batch() {
    local label="$1" root="$2"
    shift 2
    local base="${root}/Escherichia_coli"
    local sid bakta_gff prokka_gff

    echo -e "\e[32m Recording paths and staging Prokka GFFs for $# samples from ${label}... \e[0m"
    for sid in "$@"; do
        add_entry "${paths_abricate}" "Phase0_ABRicate" "${sid}" \
            "${base}/amr_Ecoli/abricate/${sid}.abricate/${sid}.abricate.card.tsv"
        add_entry "${paths_amrfinder}" "Phase0_AMRFinder" "${sid}" \
            "${base}/amr_Ecoli/amrfinder/${sid}.amrfinder/${sid}.amrfinder.tsv"
        add_entry "${paths_eggnog_bakta}" "Phase0_EggNOG_Bakta" "${sid}" \
            "${base}/annotation_Ecoli/eggnog/${sid}.eggnog/${sid}.eggnog.emapper.annotations"
        add_entry "${paths_eggnog_prokka}" "Phase0_EggNOG_Prokka" "${sid}" \
            "${base}/annotation_Ecoli/eggnog_prokka/${sid}.prokka.eggnog/${sid}.prokka.eggnog.emapper.annotations"

        bakta_gff="${base}/annotation_Ecoli/bakta/${sid}.bakta/${sid}.gff3"
        [ -s "${bakta_gff}" ] || bakta_gff="${base}/annotation_Ecoli/bakta/${sid}.bakta/${sid}.gff"
        add_entry "${paths_bakta_gff}" "Phase0_BaktaGFF" "${sid}" "${bakta_gff}"

        add_entry "${paths_quast}" "Phase0_QUAST" "${sid}" \
            "${base}/assembly_Ecoli/${sid}.filtered.quast/report.tsv"
        add_entry "${paths_isescan}" "Phase0_ISEScan" "${sid}" \
            "${base}/mge_Ecoli/isescan/${sid}.isescan"
        add_entry "${paths_integron}" "Phase0_Integron" "${sid}" \
            "${base}/mge_Ecoli/integron/${sid}.integron"
        add_entry "${paths_phispy}" "Phase0_PhiSpy" "${sid}" \
            "${base}/pai_Ecoli/phispy/${sid}.phispy"
        add_entry "${paths_ectyper}" "Phase0_ECTyper" "${sid}" \
            "${base}/typing_Ecoli/ectyper/${sid}.ectyper/output.tsv"

        # Platon plasmid classification
        platon_file="${base}/plasmid_Ecoli/platon/${sid}.platon/${sid}.tsv"
        [ -s "${platon_file}" ] || platon_file="${base}/plasmid_Ecoli/platon/${sid}.platon"
        add_entry "${paths_platon}" "Phase0_Platon" "${sid}" "${platon_file}"

        # PlasmidFinder replicon typing
        pf_file="${base}/plasmid_Ecoli/plasmidfinder/${sid}.plasmidfinder/results_tab.tsv"
        [ -s "${pf_file}" ] || pf_file="${base}/amr_Ecoli/abricate/${sid}.abricate/${sid}.abricate.plasmidfinder.tsv"
        add_entry "${paths_plasmidfinder}" "Phase0_PlasmidFinder" "${sid}" "${pf_file}"

        # MOB-suite (mob_recon contig_report.txt)
        mob_file="${base}/plasmid_Ecoli/mob_suite/${sid}.mob_recon/contig_report.txt"
        [ -s "${mob_file}" ] || mob_file="${base}/plasmid_Ecoli/mob_suite/${sid}.mob_recon"
        add_entry "${paths_mobsuite}" "Phase0_MOBsuite" "${sid}" "${mob_file}"

        # Coverage breakpoints / depth anomalies (detect_coverage_breakpoints.py)
        bp_file="${base}/assembly_Ecoli/${sid}.coverage_breakpoints/breakpoints_report.tsv"
        add_entry "${paths_breakpoints}" "Phase0_Breakpoints" "${sid}" "${bp_file}"

        # Annotation summary for ncRNAs (rRNA, tRNA, repeat_region)
        prokka_txt="${base}/annotation_Ecoli/prokka/${sid}.prokka/${sid}.prokka.txt"
        [ -s "${prokka_txt}" ] || prokka_txt="${base}/annotation_Ecoli/prokka/${sid}.prokka/${sid}.txt"
        [ -s "${prokka_txt}" ] || prokka_txt="${base}/annotation_Ecoli/bakta/${sid}.bakta/${sid}.txt"
        [ -s "${prokka_txt}" ] || prokka_txt="${prokka_gff}"
        add_entry "${paths_annotation_summary}" "Phase0_AnnotSummary" "${sid}" "${prokka_txt}"

        # Stage the Prokka GFF for Panaroo (Prokka only; Bakta .gff3 crashes Panaroo)
        prokka_gff="${base}/annotation_Ecoli/prokka/${sid}.prokka/${sid}.prokka.gff"
        [ -f "${prokka_gff}" ] || prokka_gff="${base}/annotation_Ecoli/prokka/${sid}.prokka/${sid}.gff"
        if [ ! -s "${prokka_gff}" ]; then
            echo -e "\e[31m   [${sid}] WARNING: Prokka GFF missing or empty — not staged \e[0m"
            log_failure "Phase0_StageGFF" "${sid}" "INPUT_MISSING" "Prokka GFF not found at ${prokka_gff}"
            continue
        fi
        if ! cp -f "${prokka_gff}" "${PANAROO_GFF_DIR}/${sid}.gff"; then
            echo -e "\e[31m   [${sid}] ERROR: GFF copy failed \e[0m"
            log_failure "Phase0_StageGFF" "${sid}" "EXECUTION_FAILED" "cp failed for ${prokka_gff}"
            continue
        fi
        if [ ! -s "${PANAROO_GFF_DIR}/${sid}.gff" ]; then
            log_failure "Phase0_StageGFF" "${sid}" "OUTPUT_MISSING" "Staged GFF empty: ${sid}.gff"
        fi
    done
}

# Reset the cohort list and all manifests (they are rebuilt on every run)
rm -f "${samples_list}" "${paths_abricate}" "${paths_amrfinder}" "${paths_eggnog_bakta}" \
      "${paths_eggnog_prokka}" "${paths_bakta_gff}" "${paths_quast}" "${paths_isescan}" \
      "${paths_integron}" "${paths_phispy}" "${paths_ectyper}" \
      "${paths_platon}" "${paths_plasmidfinder}" "${paths_mobsuite}" \
      "${paths_breakpoints}" "${paths_annotation_summary}"
touch "${paths_abricate}" "${paths_amrfinder}" "${paths_eggnog_bakta}" "${paths_eggnog_prokka}" \
      "${paths_bakta_gff}" "${paths_quast}" "${paths_isescan}" "${paths_integron}" \
      "${paths_phispy}" "${paths_ectyper}" \
      "${paths_platon}" "${paths_plasmidfinder}" "${paths_mobsuite}" \
      "${paths_breakpoints}" "${paths_annotation_summary}"

printf '%s\n' "${sample_ids_2[@]}" "${sample_ids_3[@]}" > "${samples_list}"

dup_ids=$(sort "${samples_list}" | uniq -d)
if [ -n "${dup_ids}" ]; then
    echo -e "\e[31m   [WARNING] Duplicate sample IDs in cohort: ${dup_ids//$'\n'/ } \e[0m"
    log_failure "Phase0_Cohort" "ALL" "DUPLICATE_ID" "Duplicate IDs: ${dup_ids//$'\n'/ }"
fi

stage_batch "bacteria_phages_2" "${BP2_PATH}" "${sample_ids_2[@]}"
stage_batch "bacteria_phages_3" "${BP3_PATH}" "${sample_ids_3[@]}"

# Automatically clean stale GFFs from earlier runs: Panaroo globs *.gff and would include them
for gff in "${PANAROO_GFF_DIR}"/*.gff; do
    [ -e "${gff}" ] || continue
    gff_id="$(basename "${gff}" .gff)"
    if ! grep -qxF "${gff_id}" "${samples_list}"; then
        echo -e "\e[33m   [STALE GFF REMOVED] ${gff_id}.gff is not in current cohort — deleted from ${PANAROO_GFF_DIR} \e[0m"
        log_failure "Phase0_StageGFF" "${gff_id}" "STALE_GFF_REMOVED" "Deleted stale GFF ${gff_id}.gff not in cohort"
        rm -f "${gff}"
    fi
done

echo -e "\e[32m Phase 0 complete: \e[0m"
echo "  Cohort list            : $(wc -l < "${samples_list}") samples in ${samples_list}"
echo "  ABRicate CARD manifest : $(wc -l < "${paths_abricate}") entries"
echo "  AMRFinderPlus manifest : $(wc -l < "${paths_amrfinder}") entries"
echo "  Bakta eggNOG manifest  : $(wc -l < "${paths_eggnog_bakta}") entries"
echo "  Prokka eggNOG manifest : $(wc -l < "${paths_eggnog_prokka}") entries"
echo "  Bakta GFF manifest     : $(wc -l < "${paths_bakta_gff}") entries"
echo "  QUAST manifest         : $(wc -l < "${paths_quast}") entries"
echo "  MGE / Plasmids / Breakpoints / Typing manifests: ${paths_isescan##*/}, ${paths_integron##*/}, ${paths_phispy##*/}, ${paths_ectyper##*/}, ${paths_platon##*/}, ${paths_mobsuite##*/}, ${paths_breakpoints##*/}, ${paths_annotation_summary##*/}"
echo "  Prokka GFFs staged     : $(find "${PANAROO_GFF_DIR}" -maxdepth 1 -name '*.gff' | wc -l) in ${PANAROO_GFF_DIR}"

########################################
# PHASE 1: AMR COLLATION (collate_amr) #
########################################

echo -e "\e[33m ====================================================== \e[0m"
echo -e "\e[33m PHASE 1: AMR COLLATION (ABRicate CARD + AMRFinderPlus)   \e[0m"
echo -e "\e[33m ====================================================== \e[0m"

echo -e "\e[31m ==================== \e[0m"
echo -e "\e[31m collate_amr.py: reconciled gene x sample matrix \e[0m"
echo -e "\e[31m ==================== \e[0m"

if [ "${force_collate}" != true ] && [ -s "${out_amr}" ]; then
    echo -e "\e[32m   [AMR] already completed — skipping \e[0m"
    log_failure "Phase1_AMR" "ALL" "SKIPPED_EXISTS" "Output ${out_amr} already exists"
elif ! manifest_has_existing "${paths_abricate}" && ! manifest_has_existing "${paths_amrfinder}"; then
    echo -e "\e[31m   [AMR] WARNING: no ABRicate or AMRFinderPlus outputs found — skipping \e[0m"
    log_failure "Phase1_AMR" "ALL" "INPUT_MISSING" "No existing file in either AMR manifest"
else
    amr_cmd=(conda run -n ncbi python3 "${collate_amr}"
        --output "${out_amr}"
        --abricate-manifest "${paths_abricate}"
        --amrfinder-manifest "${paths_amrfinder}"
        --samples "${samples_list}"
        --abricate-database "${abricate_database}"
        --per-tool-prefix "${amr_per_tool_prefix}")
    if [ "${#amrfinder_type_filter[@]}" -gt 0 ]; then
        amr_cmd+=(--amrfinder-type-filter "${amrfinder_type_filter[@]}")
    fi

    if ! "${amr_cmd[@]}"; then
        echo -e "\e[31m   [AMR] ERROR: collate_amr.py failed — continuing \e[0m"
        log_failure "Phase1_AMR" "ALL" "EXECUTION_FAILED" "collate_amr.py non-zero exit code"
    elif [ ! -s "${out_amr}" ]; then
        echo -e "\e[31m   [AMR] ERROR: expected output missing or empty after run \e[0m"
        log_failure "Phase1_AMR" "ALL" "OUTPUT_MISSING" "File ${out_amr} not created or 0 bytes"
    else
        echo -e "\e[32m   [AMR] complete -> ${out_amr} \e[0m"
    fi
fi

#######################################
# PHASE 2: PANAROO PANGENOME          #
#######################################

echo -e "\e[33m ============================================================================== \e[0m"
echo -e "\e[33m PHASE 2: PANAROO PANGENOME (ALL STAGED PROKKA GFFs) + IQTREE PHYLOGENETIC TREE    \e[0m"
echo -e "\e[33m ============================================================================== \e[0m"

panaroo_stale=false
if [ -s "${panaroo_matrix}" ]; then
    panaroo_header=$(head -n 1 "${panaroo_matrix}")
    while IFS= read -r _sid; do
        [ -n "${_sid}" ] || continue
        if ! echo "${panaroo_header}" | grep -q "${_sid}"; then
            panaroo_stale=true
            echo -e "\e[33m   [Panaroo] Cohort isolate ${_sid} missing from existing matrix — will re-run Panaroo \e[0m"
            break
        fi
    done < "${samples_list}"

    gff_count=$(find "${PANAROO_GFF_DIR}" -maxdepth 1 -name '*.gff' | wc -l)
    cohort_count=$(wc -l < "${samples_list}")
    if [ "${gff_count}" -ne "${cohort_count}" ]; then
        panaroo_stale=true
        echo -e "\e[33m   [Panaroo] Staged GFF count (${gff_count}) != cohort count (${cohort_count}) — will re-run Panaroo \e[0m"
    fi
fi

if [ "${panaroo_stale}" = true ]; then
    echo -e "\e[33m   [Panaroo] Invalidating stale pangenome outputs to ensure clean cohort analysis... \e[0m"
    rm -rf "${PANAROO_RES_DIR}"/*
    rm -f "${out_panaroo_prokka}" "${out_cluster_summary}" "${out_overview}"
fi

if [ "${panaroo_stale}" != true ] && [ -s "${panaroo_matrix}" ]; then
    echo -e "\e[32m   [Panaroo] matrix matches current cohort and already exists — skipping: ${panaroo_matrix} \e[0m"
    log_failure "Phase2_Panaroo" "ALL" "SKIPPED_EXISTS" "Output ${panaroo_matrix} already exists"
elif ! compgen -G "${PANAROO_GFF_DIR}/*.gff" > /dev/null; then
    echo -e "\e[31m   [Panaroo] ERROR: no GFF files in ${PANAROO_GFF_DIR} — skipping \e[0m"
    log_failure "Phase2_Panaroo" "ALL" "INPUT_MISSING" "No GFF files in ${PANAROO_GFF_DIR}"
else
    gff_total=$(find "${PANAROO_GFF_DIR}" -maxdepth 1 -name '*.gff' | wc -l)
    echo -e "\e[32m Running Panaroo on ${gff_total} Prokka GFFs (batch 2 + 3)... \e[0m"
    if ! conda run -n pangenome panaroo \
        -i "${PANAROO_GFF_DIR}"/*.gff \
        -o "${PANAROO_RES_DIR}" \
        --clean-mode strict \
        -a core \
        --aligner mafft \
        -t "${threads}"; then
        echo -e "\e[31m   [Panaroo] ERROR: execution failed — continuing \e[0m"
        log_failure "Phase2_Panaroo" "ALL" "EXECUTION_FAILED" "panaroo non-zero exit code"
    elif [ ! -s "${panaroo_matrix}" ]; then
        echo -e "\e[31m   [Panaroo] ERROR: expected output missing or empty after run \e[0m"
        log_failure "Phase2_Panaroo" "ALL" "OUTPUT_MISSING" "File ${panaroo_matrix} not created"
    else
        echo -e "\e[32m   [Panaroo] complete -> ${panaroo_matrix} \e[0m"
    fi
fi

# Re-evaluate alignment path in case Panaroo just created it in Phase 2
panaroo_core_aln="${PANAROO_RES_DIR}/core_gene_alignment_filtered.aln"
[ -s "${panaroo_core_aln}" ] || panaroo_core_aln="${PANAROO_RES_DIR}/core_gene_alignment.aln"

# Count samples in cohort (IQ-TREE ultrafast bootstrap -B 1000 requires >= 4 sequences)
cohort_n=$(wc -l < "${samples_list}")

if [ "${panaroo_stale:-false}" != true ] && [ "${force_collate}" != true ] && [ -s "${out_treefile}" ]; then
    echo -e "\e[32m   [IQ-TREE] treefile already exists — skipping: ${out_treefile} \e[0m"
    log_failure "Phase2.5_IQTREE" "ALL" "SKIPPED_EXISTS" "Output ${out_treefile} already exists"
elif [ ! -s "${panaroo_core_aln}" ]; then
    echo -e "\e[31m   [IQ-TREE] ERROR: Core gene alignment not found in ${PANAROO_RES_DIR} — skipping \e[0m"
    log_failure "Phase2.5_IQTREE" "ALL" "INPUT_MISSING" "No core_gene_alignment*.aln found in ${PANAROO_RES_DIR}"
elif [ "${cohort_n}" -lt 3 ]; then
    echo -e "\e[31m   [IQ-TREE] WARNING: Fewer than 3 isolates (${cohort_n}) — cannot build tree \e[0m"
    log_failure "Phase2.5_IQTREE" "ALL" "INPUT_MISSING" "Too few isolates (${cohort_n}) for tree inference"
else
    echo -e "\e[32m Running IQ-TREE on $(basename "${panaroo_core_aln}") (${cohort_n} isolates)... \e[0m"
    
    # Build IQ-TREE command (-m GTR+F+I+G4 is standard and fast for large bacterial core alignments;
    # change to -m MFP if you want full ModelFinder testing)
    iqtree_cmd=(conda run -n phylogeny iqtree
        -s "${panaroo_core_aln}"
        --prefix "${iqtree_prefix}"
        -m GTR+F+I+G4
        -T "${threads}"
        -redo)

    # Ultrafast bootstrap (-B 1000) requires at least 4 sequences
    if [ "${cohort_n}" -ge 4 ]; then
        iqtree_cmd+=(-B 1000 --alrt 1000)
    fi

    if ! "${iqtree_cmd[@]}"; then
        echo -e "\e[31m   [IQ-TREE] ERROR: execution failed — continuing \e[0m"
        log_failure "Phase2.5_IQTREE" "ALL" "EXECUTION_FAILED" "iqtree non-zero exit code"
    elif [ ! -s "${out_treefile}" ]; then
        echo -e "\e[31m   [IQ-TREE] ERROR: expected .treefile missing or empty after run \e[0m"
        log_failure "Phase2.5_IQTREE" "ALL" "OUTPUT_MISSING" "File ${out_treefile} not created"
    else
        echo -e "\e[32m   [IQ-TREE] complete -> ${out_treefile} \e[0m"
    fi
fi

#################################################
# PHASE 3: ANNOTATION COLLATION (collate_annot) #
#################################################

echo -e "\e[33m ====================================================== \e[0m"
echo -e "\e[33m PHASE 3: ANNOTATION COLLATION (GFF + Panaroo + eggNOG)   \e[0m"
echo -e "\e[33m ====================================================== \e[0m"

# Step 3.1: Bakta GFF3 + Bakta eggNOG (GFF mode, no pangenome)
echo -e "\e[31m ==================== \e[0m"
echo -e "\e[31m collate_annotations.py [3.1]: Bakta GFF + Bakta eggNOG \e[0m"
echo -e "\e[31m ==================== \e[0m"

if [ "${force_collate}" != true ] && [ -s "${out_bakta_eggnog}" ]; then
    echo -e "\e[32m   [Bakta+eggNOG] already completed — skipping \e[0m"
    log_failure "Phase3_Bakta_EggNOG" "ALL" "SKIPPED_EXISTS" "Output ${out_bakta_eggnog} exists"
elif ! manifest_has_existing "${paths_bakta_gff}" || ! manifest_has_existing "${paths_eggnog_bakta}"; then
    echo -e "\e[31m   [Bakta+eggNOG] WARNING: Bakta GFFs or Bakta eggNOG files missing — skipping \e[0m"
    log_failure "Phase3_Bakta_EggNOG" "ALL" "INPUT_MISSING" "No existing Bakta GFF or Bakta eggNOG file"
else
    if ! conda run -n ncbi python3 "${collate_annotations}" \
        --output "${out_bakta_eggnog}" \
        --gff-manifest "${paths_bakta_gff}" \
        --eggnog-manifest "${paths_eggnog_bakta}" \
        --samples "${samples_list}" \
        --eggnog-summary cog kegg completeness \
        --eggnog-summary-prefix "${eggnog_bakta_prefix}"; then
        echo -e "\e[31m   [Bakta+eggNOG] ERROR: collate_annotations.py failed — continuing \e[0m"
        log_failure "Phase3_Bakta_EggNOG" "ALL" "EXECUTION_FAILED" "collate_annotations.py non-zero exit"
    elif [ ! -s "${out_bakta_eggnog}" ]; then
        echo -e "\e[31m   [Bakta+eggNOG] ERROR: expected output missing or empty after run \e[0m"
        log_failure "Phase3_Bakta_EggNOG" "ALL" "OUTPUT_MISSING" "File ${out_bakta_eggnog} not created"
    else
        echo -e "\e[32m   [Bakta+eggNOG] complete -> ${out_bakta_eggnog} \e[0m"
    fi
fi

# Step 3.2: Panaroo clusters + Prokka GFF (+ Prokka eggNOG when available)
echo -e "\e[31m ==================== \e[0m"
echo -e "\e[31m collate_annotations.py [3.2]: Panaroo + Prokka (+ eggNOG) \e[0m"
echo -e "\e[31m ==================== \e[0m"

if [ "${force_collate}" != true ] && [ -s "${out_panaroo_prokka}" ]; then
    echo -e "\e[32m   [Panaroo+Prokka] already completed — skipping \e[0m"
    log_failure "Phase3_Panaroo_Prokka" "ALL" "SKIPPED_EXISTS" "Output ${out_panaroo_prokka} exists"
elif [ ! -s "${panaroo_matrix}" ]; then
    echo -e "\e[31m   [Panaroo+Prokka] WARNING: Panaroo matrix missing — skipping \e[0m"
    log_failure "Phase3_Panaroo_Prokka" "ALL" "INPUT_MISSING" "Panaroo matrix not found at ${panaroo_matrix}"
else
    panaroo_annot_cmd=(conda run -n ncbi python3 "${collate_annotations}"
        --output "${out_panaroo_prokka}"
        --panaroo "${panaroo_matrix}"
        --gff-dir "${PANAROO_GFF_DIR}"
        --samples "${samples_list}"
        --cluster-summary "${out_cluster_summary}"
        --core-cutoff "${core_cutoff}"
        --soft-core-cutoff "${soft_core_cutoff}"
        --shell-cutoff "${shell_cutoff}")

    # Attach Prokka eggNOG only if at least one file exists (eggNOG on Prokka = pending otherwise)
    if manifest_has_existing "${paths_eggnog_prokka}"; then
        panaroo_annot_cmd+=(--eggnog-manifest "${paths_eggnog_prokka}"
            --eggnog-summary cog kegg completeness
            --eggnog-summary-prefix "${eggnog_prokka_prefix}")
    else
        echo -e "\e[33m   [Panaroo+Prokka] Prokka eggNOG files not available — running without eggNOG \e[0m"
        log_failure "Phase3_Panaroo_Prokka" "ALL" "INPUT_MISSING" "No Prokka eggNOG files; eggNOG omitted"
    fi

    if ! "${panaroo_annot_cmd[@]}"; then
        echo -e "\e[31m   [Panaroo+Prokka] ERROR: collate_annotations.py failed — continuing \e[0m"
        log_failure "Phase3_Panaroo_Prokka" "ALL" "EXECUTION_FAILED" "collate_annotations.py non-zero exit"
    elif [ ! -s "${out_panaroo_prokka}" ]; then
        echo -e "\e[31m   [Panaroo+Prokka] ERROR: expected output missing or empty after run \e[0m"
        log_failure "Phase3_Panaroo_Prokka" "ALL" "OUTPUT_MISSING" "File ${out_panaroo_prokka} not created"
    else
        echo -e "\e[32m   [Panaroo+Prokka] complete -> ${out_panaroo_prokka} \e[0m"
    fi
fi

#################################################
# PHASE 4: MASTER GENOME OVERVIEW (QC + typing) #
#################################################

echo -e "\e[33m ====================================================== \e[0m"
echo -e "\e[33m PHASE 4: MASTER GENOME OVERVIEW (QC + TYPING + MGE)      \e[0m"
echo -e "\e[33m ====================================================== \e[0m"

echo -e "\e[31m ==================== \e[0m"
echo -e "\e[31m collate_genome_overview.py: QUAST+CheckM+Panaroo+MLST+ECTyper+MGE \e[0m"
echo -e "\e[31m ==================== \e[0m"

if [ "${force_collate}" != true ] && [ -s "${out_overview}" ]; then
    echo -e "\e[32m   [Overview] already completed — skipping \e[0m"
    log_failure "Phase4_Overview" "ALL" "SKIPPED_EXISTS" "Output ${out_overview} already exists"
else
    # Batch-wide CheckM / MLST tables: pass whichever exist
    checkm_inputs=()
    [ -s "${checkm_bp2}" ] && checkm_inputs+=("${checkm_bp2}")
    [ -s "${checkm_bp3}" ] && checkm_inputs+=("${checkm_bp3}")
    mlst_inputs=()
    [ -s "${mlst_bp2}" ] && mlst_inputs+=("${mlst_bp2}")
    [ -s "${mlst_bp3}" ] && mlst_inputs+=("${mlst_bp3}")

    overview_cmd=(conda run -n ncbi python3 "${collate_genome_overview}"
        --output "${out_overview}"
        --samples "${samples_list}"
        --quast-manifest "${paths_quast}"
        --ectyper "${paths_ectyper}"
        --isescan-manifest "${paths_isescan}"
        --integron-manifest "${paths_integron}"
        --phispy-manifest "${paths_phispy}"
        --detail-prefix "${overview_detail_prefix}"
        --species "${species}"
        --size-status
        --min-completeness "${min_completeness}"
        --max-contamination "${max_contamination}"
        --max-contigs "${max_contigs}"
        --core-cutoff "${core_cutoff}"
        --soft-core-cutoff "${soft_core_cutoff}"
        --shell-cutoff "${shell_cutoff}")

    [ "${#checkm_inputs[@]}" -gt 0 ] && overview_cmd+=(--checkm "${checkm_inputs[@]}")
    [ "${#mlst_inputs[@]}" -gt 0 ] && overview_cmd+=(--mlst "${mlst_inputs[@]}")
    if [ -s "${panaroo_matrix}" ]; then
        overview_cmd+=(--panaroo "${panaroo_matrix}")
    else
        echo -e "\e[33m   [Overview] Panaroo matrix not found — pangenome columns will be empty \e[0m"
        log_failure "Phase4_Overview" "ALL" "INPUT_MISSING" "Panaroo matrix not found at ${panaroo_matrix}"
    fi

    # Conditionally attach Plasmid, Breakpoint, and Annotation manifests if outputs exist
    if manifest_has_existing "${paths_platon}"; then
        overview_cmd+=(--platon-manifest "${paths_platon}")
    fi
    if manifest_has_existing "${paths_plasmidfinder}"; then
        overview_cmd+=(--plasmidfinder-manifest "${paths_plasmidfinder}")
    fi
    if manifest_has_existing "${paths_mobsuite}"; then
        overview_cmd+=(--mobsuite-manifest "${paths_mobsuite}")
    fi
    if manifest_has_existing "${paths_breakpoints}"; then
        overview_cmd+=(--breakpoints-manifest "${paths_breakpoints}")
    fi
    if manifest_has_existing "${paths_annotation_summary}"; then
        overview_cmd+=(--annotation-manifest "${paths_annotation_summary}")
    fi

    [ "${#checkm_inputs[@]}" -eq 0 ] && \
        log_failure "Phase4_Overview" "ALL" "INPUT_MISSING" "No CheckM table in batch 2 or 3"
    [ "${#mlst_inputs[@]}" -eq 0 ] && \
        log_failure "Phase4_Overview" "ALL" "INPUT_MISSING" "No MLST table in batch 2 or 3"

    if ! "${overview_cmd[@]}"; then
        echo -e "\e[31m   [Overview] ERROR: collate_genome_overview.py failed — continuing \e[0m"
        log_failure "Phase4_Overview" "ALL" "EXECUTION_FAILED" "collate_genome_overview.py non-zero exit"
    elif [ ! -s "${out_overview}" ]; then
        echo -e "\e[31m   [Overview] ERROR: expected output missing or empty after run \e[0m"
        log_failure "Phase4_Overview" "ALL" "OUTPUT_MISSING" "File ${out_overview} not created"
    else
        echo -e "\e[32m   [Overview] complete -> ${out_overview} \e[0m"
    fi
fi

###########################################################
# SUMMARY REPORT                                          #
###########################################################

fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING\]\|\[INPUT_MISSING\]' \
    "${FAIL_LOG}" 2>/dev/null || true)
skip_count=$(grep -c '\[SKIPPED' "${FAIL_LOG}" 2>/dev/null || true)

echo ""
echo -e "\e[32m ======================================================== \e[0m"
echo -e "\e[32m E. COLI TN MERGE COMPLETE — $(date) \e[0m"
echo -e "\e[32m ======================================================== \e[0m"
echo ""
echo -e "\e[32m -- COHORT & MANIFESTS ------------------------------------ \e[0m"
echo -e "\e[32m  Cohort list        : ${samples_list} \e[0m"
echo -e "\e[32m  Manifests          : ${PATHS_DIR} \e[0m"
echo ""
echo -e "\e[32m -- PHASE 1: AMR ------------------------------------------ \e[0m"
echo -e "\e[32m  Reconciled matrix  : ${out_amr} \e[0m"
echo -e "\e[32m  Per-tool matrices  : ${amr_per_tool_prefix}_{abricate,amrfinder}_matrix.tsv \e[0m"
echo ""
echo -e "\e[32m -- PHASE 2: PANGENOME ------------------------------------ \e[0m"
echo -e "\e[32m  Panaroo matrix     : ${panaroo_matrix} \e[0m"
echo ""
echo -e "\e[32m -- PHASE 3: ANNOTATION ----------------------------------- \e[0m"
echo -e "\e[32m  Bakta + eggNOG     : ${out_bakta_eggnog} \e[0m"
echo -e "\e[32m  Panaroo + Prokka   : ${out_panaroo_prokka} \e[0m"
echo -e "\e[32m  Cluster summary    : ${out_cluster_summary} \e[0m"
echo -e "\e[32m  eggNOG summaries   : ${INDIVIDUAL_OUT_DIR}/eggnog_{bakta,prokka}_{cog,kegg,completeness}.tsv \e[0m"
echo ""
echo -e "\e[32m -- PHASE 4: GENOME OVERVIEW ------------------------------ \e[0m"
echo -e "\e[32m  Overview QC table  : ${out_overview} \e[0m"
echo -e "\e[32m  Detail tables      : ${overview_detail_prefix}_*.tsv \e[0m"
echo ""
echo -e "\e[32m -- LOG FILES --------------------------------------------- \e[0m"
echo -e "\e[32m  Full execution log : ${LOG} \e[0m"

if [ "${fail_count}" -gt 0 ]; then
    echo -e "\e[31m  Failure/Issues log : ${FAIL_LOG} (${fail_count} failures/missing inputs!) \e[0m"
    echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which samples/steps failed and why. \e[0m"
else
    echo -e "\e[32m  Failure/Issues log : ${FAIL_LOG} (0 errors recorded) \e[0m"
fi
echo -e "\e[32m  Steps skipped (outputs already present): ${skip_count} \e[0m"

echo ""
echo -e "\e[32m -- SUGGESTED NEXT STEPS ---------------------------------- \e[0m"
echo -e "\e[32m  1. Read the *.summary.txt next to each collated table (dropped/unmatched isolates) \e[0m"
echo -e "\e[32m  2. Check the eggNOG match rate in the annotation summaries (Prokka vs Bakta mismatch) \e[0m"
echo -e "\e[32m  3. Re-run with force_collate=true after changing the cohort, thresholds or inputs \e[0m"
echo -e "\e[32m  Full log: ${LOG} \e[0m"

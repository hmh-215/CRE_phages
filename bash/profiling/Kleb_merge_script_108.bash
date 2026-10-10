#!/bin/bash
# Merge all Hospital 108 K. pneumoniae datasets across batch 1, batch 2 and batch 3
# Workflow (collate_*.py edition):
#   Phase 0 : Cohort list + per-sample manifests (sample<TAB>path) + stage Prokka .gff files
#             (Prokka only: Panaroo crashes on Bakta .gff3)
#   Phase 1 : AMR collation        -> collate_amr.py (ABRicate CARD + AMRFinderPlus)
#   Phase 2 : Panaroo pangenome    -> panaroo_Kleb_108/results/
#   Phase 3 : Annotation collation -> collate_annotations.py
#             - Step 3.1: Bakta GFF3 + Bakta eggNOG        (GFF mode)
#             - Step 3.2: Panaroo + Prokka GFF + Prokka eggNOG (Panaroo mode)
#   Phase 4 : Genome overview      -> collate_genome_overview.py
#             (QUAST + CheckM + Panaroo + MLST + Kleborate + ISEScan/IntegronFinder/PhiSpy)
# Building No. 2
set -euo pipefail

################
# GLOBAL SETUP #
################

# Base project locations (with fallback handling for server path variations)
BP1_PATH="/storage/student9/projects/bacteria_phages"
[ ! -d "${BP1_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage" ] && \
    BP1_PATH="/storage/student9/projects/bacteria_phage"
[ ! -d "${BP1_PATH}" ] && [ -d "/storage/student9/bacteria_phages" ] && \
    BP1_PATH="/storage/student9/bacteria_phages"

BP2_PATH="/storage/student9/projects/bacteria_phages_2"
[ ! -d "${BP2_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage_2" ] && \
    BP2_PATH="/storage/student9/projects/bacteria_phage_2"
[ ! -d "${BP2_PATH}" ] && [ -d "/storage/student9/bacteria_phages_2" ] && \
    BP2_PATH="/storage/student9/bacteria_phages_2"

BP3_PATH="/storage/student9/projects/bacteria_phages_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage_3" ] && \
    BP3_PATH="/storage/student9/projects/bacteria_phage_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/bacteria_phages_3" ] && \
    BP3_PATH="/storage/student9/bacteria_phages_3"
[ ! -d "${BP3_PATH}" ] && [ -d "/storage/student9/bacteria_phage_3" ] && \
    BP3_PATH="/storage/student9/bacteria_phage_3"

# Target reorganization root folder in bacteria_phages_3
MERGE_PATH="${BP3_PATH}/merge_Kpneumoniae_108"

# Subfolder organization:
# merge_Kpneumoniae_108/
#   |___ paths/                    (cohort list + manifests)
#   |___ panaroo_Kleb_108/
#   |     |___ gff/                (staged Prokka GFFs)
#   |     |___ results/            (Panaroo outputs)
#   |___ merge_individual_outputs/ (per-tool matrices, eggNOG summaries, overview details)
#   |___ merge_metadata/           (final collated tables)
PATHS_DIR="${MERGE_PATH}/paths"
PANAROO_DIR="${MERGE_PATH}/panaroo_Kleb_108"
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

# Cohort: Hospital 108 isolates across batch 1, 2, and 3
sample_ids_1=(
    WS2762512A01 WS2762512A02 WS2762512A03 WS2762512A04 WS2762512A06
    WS2762512A07 WS2762512A10 WS2762512A14 WS2762512A15
)
sample_ids_2=(
    K10-LJF30359 K14-LJF30358 K16-LJF30357 K66-LJF30363 K74-LJF30355 K75-LJF30356
)
sample_ids_3=(
    WS2762607A54
)

# Cohort list + manifests (Phase 0 writes them; every collate_*.py reads them)
# Manifest format: one "sample<TAB>path" per line; the cohort file is one sample per line.
samples_list="${PATHS_DIR}/samples_108.txt"
paths_abricate="${PATHS_DIR}/paths_abricate_card_108.txt"
paths_amrfinder="${PATHS_DIR}/paths_amrfinder_108.txt"
paths_eggnog_bakta="${PATHS_DIR}/paths_eggnog_bakta_108.txt"
paths_eggnog_prokka="${PATHS_DIR}/paths_eggnog_prokka_108.txt"
paths_bakta_gff="${PATHS_DIR}/paths_bakta_gff_108.txt"
paths_quast="${PATHS_DIR}/paths_quast_108.txt"
paths_isescan="${PATHS_DIR}/paths_isescan_108.txt"
paths_integron="${PATHS_DIR}/paths_integron_108.txt"
paths_phispy="${PATHS_DIR}/paths_phispy_108.txt"
paths_platon="${PATHS_DIR}/paths_platon_108.txt"
paths_plasmidfinder="${PATHS_DIR}/paths_plasmidfinder_108.txt"
paths_mobsuite="${PATHS_DIR}/paths_mobsuite_108.txt"
paths_breakpoints="${PATHS_DIR}/paths_breakpoints_108.txt"
paths_annotation_summary="${PATHS_DIR}/paths_annotation_summary_108.txt"

# Batch-wide tables passed straight to collate_genome_overview.py (no manifest needed)
checkm_bp1="${BP1_PATH}/Klebsiella_pneumoniae/annotation_Kleb/checkm/all_Kleb.quality.checkm.tsv"
[ -s "${checkm_bp1}" ] || checkm_bp1="${BP1_PATH}/Klebsiella_pneumoniae/reassembly_Kleb/checkm/Kleb.reassembly.quality.checkm.tsv"

checkm_bp2="${BP2_PATH}/Klebsiella_pneumoniae/annotation_Kleb/checkm/all_Kleb.quality.checkm.tsv"
[ -s "${checkm_bp2}" ] || checkm_bp2="${BP2_PATH}/Klebsiella_pneumoniae/reassembly_Kleb/checkm/Kleb.reassembly.quality.checkm.tsv"

checkm_bp3="${BP3_PATH}/Klebsiella_pneumoniae/annotation_Kleb/checkm/all_Kleb.quality.checkm.tsv"
[ -s "${checkm_bp3}" ] || checkm_bp3="${BP3_PATH}/Klebsiella_pneumoniae/reassembly_Kleb/checkm/Kleb.reassembly.quality.checkm.tsv"

mlst_bp1="${BP1_PATH}/Klebsiella_pneumoniae/annotation_Kleb/mlst/all_Kleb.mlst.csv"
[ -s "${mlst_bp1}" ] || mlst_bp1="${BP1_PATH}/Klebsiella_pneumoniae/annotation_Kleb/mlst/all_Kleb.mlst.tsv"

mlst_bp2="${BP2_PATH}/Klebsiella_pneumoniae/annotation_Kleb/mlst/all_Kleb.mlst.csv"
[ -s "${mlst_bp2}" ] || mlst_bp2="${BP2_PATH}/Klebsiella_pneumoniae/annotation_Kleb/mlst/all_Kleb.mlst.tsv"

mlst_bp3="${BP3_PATH}/Klebsiella_pneumoniae/annotation_Kleb/mlst/all_Kleb.mlst.csv"
[ -s "${mlst_bp3}" ] || mlst_bp3="${BP3_PATH}/Klebsiella_pneumoniae/annotation_Kleb/mlst/all_Kleb.mlst.tsv"

kleborate_bp1="${BP1_PATH}/Klebsiella_pneumoniae/pai_Kleb/kleborate/klebsiella_pneumo_complex_output.txt"
[ -s "${kleborate_bp1}" ] || kleborate_bp1="${BP1_PATH}/Klebsiella_pneumoniae/pai_Kleb/kleborate/all_Kleb.kleborate.tsv"

kleborate_bp2="${BP2_PATH}/Klebsiella_pneumoniae/pai_Kleb/kleborate/klebsiella_pneumo_complex_output.txt"
[ -s "${kleborate_bp2}" ] || kleborate_bp2="${BP2_PATH}/Klebsiella_pneumoniae/pai_Kleb/kleborate/all_Kleb.kleborate.tsv"

kleborate_bp3="${BP3_PATH}/Klebsiella_pneumoniae/pai_Kleb/kleborate/klebsiella_pneumo_complex_output.txt"
[ -s "${kleborate_bp3}" ] || kleborate_bp3="${BP3_PATH}/Klebsiella_pneumoniae/pai_Kleb/kleborate/all_Kleb.kleborate.tsv"

# Output files
panaroo_matrix="${PANAROO_RES_DIR}/gene_presence_absence.csv"
out_amr="${METADATA_OUT_DIR}/amr_reconciled_108.tsv"
amr_per_tool_prefix="${INDIVIDUAL_OUT_DIR}/amr_108"
out_bakta_eggnog="${METADATA_OUT_DIR}/bakta_eggnog_annotated_genes_108.tsv"
eggnog_bakta_prefix="${INDIVIDUAL_OUT_DIR}/eggnog_bakta"
out_panaroo_prokka="${METADATA_OUT_DIR}/panaroo_prokka_annotated_genes_108.tsv"
out_cluster_summary="${METADATA_OUT_DIR}/panaroo_prokka_clusters_summary_108.tsv"
eggnog_prokka_prefix="${INDIVIDUAL_OUT_DIR}/eggnog_prokka"
out_overview="${METADATA_OUT_DIR}/genome_overview_qc_108.tsv"
overview_detail_prefix="${INDIVIDUAL_OUT_DIR}/overview_108"

# Tunable parameters
threads=16                  # Panaroo threads
force_collate=false         # true = regenerate collated tables even if they already exist
abricate_database="card"    # ABRicate database label
amrfinder_type_filter=()    # e.g. (AMR) keeps only Type==AMR rows; empty = keep all types
core_cutoff=0.99            # Panaroo compartments: Core    >= this cohort frequency
soft_core_cutoff=0.95       #                       SoftCore >= this
shell_cutoff=0.15           #                       Shell    >= this (below = Cloud)
species="klebsiella"        # overview size range: ecoli | klebsiella | pseudomonas
min_completeness=95.0       # CheckM PASS threshold (%)
max_contamination=5.0       # CheckM PASS threshold (%)
max_contigs=300             # Assembly_Status PASS threshold

# Custom Python helper scripts (priority: bacteria_phage_3/scripts, fallback: /storage/student9/tools)
TOOL_PATH="${BP3_PATH}/scripts"
[ ! -d "${TOOL_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phages_3/scripts" ] && \
    TOOL_PATH="/storage/student9/projects/bacteria_phages_3/scripts"
[ ! -d "${TOOL_PATH}" ] && [ -d "/storage/student9/projects/bacteria_phage_3/scripts" ] && \
    TOOL_PATH="/storage/student9/projects/bacteria_phage_3/scripts"
[ ! -d "${TOOL_PATH}" ] && [ -d "/storage/student9/tools" ] && \
    TOOL_PATH="/storage/student9/tools"

collate_amr="${TOOL_PATH}/collate_amr.py"
collate_annotations="${TOOL_PATH}/collate_annotations.py"
collate_genome_overview="${TOOL_PATH}/collate_genome_overview.py"

# Dual-log setup
LOG="${MERGE_PATH}/merge_script_108.log"
FAIL_LOG="${MERGE_PATH}/merge_script_108.failed_skipped.log"
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
echo " Merge all Hospital 108 K. pneumoniae samples — Started: $(date)"
echo " Cohort         : ${#sample_ids_1[@]} (batch 1) + ${#sample_ids_2[@]} (batch 2) + ${#sample_ids_3[@]} (batch 3) isolates"
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
stage_batch() {
    local label="$1" root="$2"
    shift 2
    local base="${root}/Klebsiella_pneumoniae"
    local sid bakta_gff prokka_gff quast_rpt

    echo -e "\e[32m Recording paths and staging Prokka GFFs for $# samples from ${label}... \e[0m"
    for sid in "$@"; do
        add_entry "${paths_abricate}" "Phase0_ABRicate" "${sid}" \
            "${base}/amr_Kleb/abricate/${sid}.abricate/${sid}.abricate.card.tsv"
        add_entry "${paths_amrfinder}" "Phase0_AMRFinder" "${sid}" \
            "${base}/annotation_Kleb/amrfinder/${sid}.amrfinder/${sid}.amrfinder.tsv"
        add_entry "${paths_eggnog_bakta}" "Phase0_EggNOG_Bakta" "${sid}" \
            "${base}/annotation_Kleb/eggnog/${sid}.eggnog/${sid}.eggnog.emapper.annotations"
        add_entry "${paths_eggnog_prokka}" "Phase0_EggNOG_Prokka" "${sid}" \
            "${base}/annotation_Kleb/eggnog_prokka/${sid}.prokka.eggnog/${sid}.prokka.eggnog.emapper.annotations"

        bakta_gff="${base}/annotation_Kleb/bakta/${sid}.bakta/${sid}.bakta.gff3"
        [ -s "${bakta_gff}" ] || bakta_gff="${base}/annotation_Kleb/bakta/${sid}.bakta/${sid}.gff3"
        [ -s "${bakta_gff}" ] || bakta_gff="${base}/annotation_Kleb/bakta/${sid}.bakta/${sid}.bakta.gff"
        [ -s "${bakta_gff}" ] || bakta_gff="${base}/annotation_Kleb/bakta/${sid}.bakta/${sid}.gff"
        add_entry "${paths_bakta_gff}" "Phase0_BaktaGFF" "${sid}" "${bakta_gff}"

        quast_rpt="${base}/assembly_Kleb/${sid}.filtered.quast/report.tsv"
        [ -s "${quast_rpt}" ] || quast_rpt="${base}/reassembly_Kleb/${sid}.unicycler/quast/report.tsv"
        add_entry "${paths_quast}" "Phase0_QUAST" "${sid}" "${quast_rpt}"

        add_entry "${paths_isescan}" "Phase0_ISEScan" "${sid}" \
            "${base}/mge_Kleb/isescan/${sid}.isescan"
        add_entry "${paths_integron}" "Phase0_Integron" "${sid}" \
            "${base}/mge_Kleb/integron/${sid}.integron"
        add_entry "${paths_phispy}" "Phase0_PhiSpy" "${sid}" \
            "${base}/pai_Kleb/phispy/${sid}.phispy"

        # Platon plasmid classification
        platon_file="${base}/plasmid_Kleb/platon/${sid}.platon/${sid}.tsv"
        [ -s "${platon_file}" ] || platon_file="${base}/plasmid_Kleb/platon/${sid}.platon"
        add_entry "${paths_platon}" "Phase0_Platon" "${sid}" "${platon_file}"

        # PlasmidFinder replicon typing
        pf_file="${base}/plasmid_Kleb/plasmidfinder/${sid}.plasmidfinder/results_tab.tsv"
        [ -s "${pf_file}" ] || pf_file="${base}/amr_Kleb/abricate/${sid}.abricate/${sid}.abricate.plasmidfinder.tsv"
        add_entry "${paths_plasmidfinder}" "Phase0_PlasmidFinder" "${sid}" "${pf_file}"

        # MOB-suite (mob_recon contig_report.txt)
        mob_file="${base}/plasmid_Kleb/mob_suite/${sid}.mob_recon/contig_report.txt"
        [ -s "${mob_file}" ] || mob_file="${base}/plasmid_Kleb/mob_suite/${sid}.mob_recon"
        add_entry "${paths_mobsuite}" "Phase0_MOBsuite" "${sid}" "${mob_file}"

        # Coverage breakpoints / depth anomalies (detect_coverage_breakpoints.py)
        bp_file="${base}/assembly_Kleb/${sid}.coverage_breakpoints/breakpoints_report.tsv"
        add_entry "${paths_breakpoints}" "Phase0_Breakpoints" "${sid}" "${bp_file}"

        # Annotation summary for ncRNAs (rRNA, tRNA, repeat_region)
        prokka_txt="${base}/annotation_Kleb/prokka/${sid}.prokka/${sid}.prokka.txt"
        [ -s "${prokka_txt}" ] || prokka_txt="${base}/annotation_Kleb/prokka/${sid}.prokka/${sid}.txt"
        [ -s "${prokka_txt}" ] || prokka_txt="${base}/annotation_Kleb/bakta/${sid}.bakta/${sid}.txt"
        [ -s "${prokka_txt}" ] || prokka_txt="${prokka_gff}"
        add_entry "${paths_annotation_summary}" "Phase0_AnnotSummary" "${sid}" "${prokka_txt}"

        # Stage the Prokka GFF for Panaroo (Prokka only; Bakta .gff3 crashes Panaroo)
        prokka_gff="${base}/annotation_Kleb/prokka/${sid}.prokka/${sid}.prokka.gff"
        [ -f "${prokka_gff}" ] || prokka_gff="${base}/annotation_Kleb/prokka/${sid}.prokka/${sid}.gff"
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

# Reset the cohort list and all manifests
rm -f "${samples_list}" "${paths_abricate}" "${paths_amrfinder}" "${paths_eggnog_bakta}" \
      "${paths_eggnog_prokka}" "${paths_bakta_gff}" "${paths_quast}" "${paths_isescan}" \
      "${paths_integron}" "${paths_phispy}" \
      "${paths_platon}" "${paths_plasmidfinder}" "${paths_mobsuite}" \
      "${paths_breakpoints}" "${paths_annotation_summary}"
touch "${paths_abricate}" "${paths_amrfinder}" "${paths_eggnog_bakta}" "${paths_eggnog_prokka}" \
      "${paths_bakta_gff}" "${paths_quast}" "${paths_isescan}" "${paths_integron}" \
      "${paths_phispy}" \
      "${paths_platon}" "${paths_plasmidfinder}" "${paths_mobsuite}" \
      "${paths_breakpoints}" "${paths_annotation_summary}"

printf '%s\n' "${sample_ids_1[@]}" "${sample_ids_2[@]}" "${sample_ids_3[@]}" > "${samples_list}"

dup_ids=$(sort "${samples_list}" | uniq -d)
if [ -n "${dup_ids}" ]; then
    echo -e "\e[31m   [WARNING] Duplicate sample IDs in cohort: ${dup_ids//$'\n'/ } \e[0m"
    log_failure "Phase0_Cohort" "ALL" "DUPLICATE_ID" "Duplicate IDs: ${dup_ids//$'\n'/ }"
fi

stage_batch "bacteria_phages (batch 1)" "${BP1_PATH}" "${sample_ids_1[@]}"
stage_batch "bacteria_phages_2 (batch 2)" "${BP2_PATH}" "${sample_ids_2[@]}"
stage_batch "bacteria_phages_3 (batch 3)" "${BP3_PATH}" "${sample_ids_3[@]}"

# Warn about stale GFFs from earlier runs
for gff in "${PANAROO_GFF_DIR}"/*.gff; do
    [ -e "${gff}" ] || continue
    gff_id="$(basename "${gff}" .gff)"
    if ! grep -qxF "${gff_id}" "${samples_list}"; then
        echo -e "\e[31m   [WARNING] ${gff_id}.gff is not in the cohort — remove it before Panaroo \e[0m"
        log_failure "Phase0_StageGFF" "${gff_id}" "STALE_GFF" "GFF in ${PANAROO_GFF_DIR} not in cohort"
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
echo "  MGE / Plasmids / Breakpoints / Typing manifests: ${paths_isescan##*/}, ${paths_integron##*/}, ${paths_phispy##*/}, ${paths_platon##*/}, ${paths_mobsuite##*/}, ${paths_breakpoints##*/}, ${paths_annotation_summary##*/}"
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

echo -e "\e[33m ====================================================== \e[0m"
echo -e "\e[33m PHASE 2: PANAROO PANGENOME (ALL STAGED PROKKA GFFs)     \e[0m"
echo -e "\e[33m ====================================================== \e[0m"

if [ -s "${panaroo_matrix}" ]; then
    echo -e "\e[32m   [Panaroo] matrix already exists — skipping: ${panaroo_matrix} \e[0m"
    log_failure "Phase2_Panaroo" "ALL" "SKIPPED_EXISTS" "Output ${panaroo_matrix} already exists"
elif ! compgen -G "${PANAROO_GFF_DIR}/*.gff" > /dev/null; then
    echo -e "\e[31m   [Panaroo] ERROR: no GFF files in ${PANAROO_GFF_DIR} — skipping \e[0m"
    log_failure "Phase2_Panaroo" "ALL" "INPUT_MISSING" "No GFF files in ${PANAROO_GFF_DIR}"
else
    gff_total=$(find "${PANAROO_GFF_DIR}" -maxdepth 1 -name '*.gff' | wc -l)
    echo -e "\e[32m Running Panaroo on ${gff_total} Prokka GFFs across batches 1, 2, and 3... \e[0m"
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
echo -e "\e[31m collate_genome_overview.py: QUAST+CheckM+Panaroo+MLST+Kleborate+MGE \e[0m"
echo -e "\e[31m ==================== \e[0m"

if [ "${force_collate}" != true ] && [ -s "${out_overview}" ]; then
    echo -e "\e[32m   [Overview] already completed — skipping \e[0m"
    log_failure "Phase4_Overview" "ALL" "SKIPPED_EXISTS" "Output ${out_overview} already exists"
else
    # Batch-wide CheckM tables: pass whichever exist
    checkm_inputs=()
    [ -s "${checkm_bp1}" ] && checkm_inputs+=("${checkm_bp1}")
    [ -s "${checkm_bp2}" ] && checkm_inputs+=("${checkm_bp2}")
    [ -s "${checkm_bp3}" ] && checkm_inputs+=("${checkm_bp3}")

    # MLST tables: pass whichever exist
    mlst_inputs=()
    [ -s "${mlst_bp1}" ] && mlst_inputs+=("${mlst_bp1}")
    [ -s "${mlst_bp2}" ] && mlst_inputs+=("${mlst_bp2}")
    [ -s "${mlst_bp3}" ] && mlst_inputs+=("${mlst_bp3}")

    # Kleborate tables: pass whichever exist
    kleborate_inputs=()
    [ -s "${kleborate_bp1}" ] && kleborate_inputs+=("${kleborate_bp1}")
    [ -s "${kleborate_bp2}" ] && kleborate_inputs+=("${kleborate_bp2}")
    [ -s "${kleborate_bp3}" ] && kleborate_inputs+=("${kleborate_bp3}")

    overview_cmd=(conda run -n ncbi python3 "${collate_genome_overview}"
        --output "${out_overview}"
        --samples "${samples_list}"
        --quast-manifest "${paths_quast}"
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
    [ "${#kleborate_inputs[@]}" -gt 0 ] && overview_cmd+=(--kleborate "${kleborate_inputs[@]}")

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
        log_failure "Phase4_Overview" "ALL" "INPUT_MISSING" "No CheckM table in batches 1, 2, or 3"
    [ "${#kleborate_inputs[@]}" -eq 0 ] && [ "${#mlst_inputs[@]}" -eq 0 ] && \
        log_failure "Phase4_Overview" "ALL" "INPUT_MISSING" "No Kleborate or MLST table in batches 1, 2, or 3"

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
echo -e "\e[32m K. PNEUMONIAE HOSPITAL 108 MERGE COMPLETE — $(date)      \e[0m"
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
echo -e "\e[32m  Staged Prokka GFFs : ${PANAROO_GFF_DIR}/ \e[0m"
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
echo -e "\e[32m  1. Read *.summary.txt next to each collated table (dropped/unmatched isolates) \e[0m"
echo -e "\e[32m  2. Inspect Panaroo gene_presence_absence.csv for core genome estimation \e[0m"
echo -e "\e[32m  3. Re-run with force_collate=true if inputs, cohorts, or thresholds change \e[0m"
echo -e "\e[32m  Full log: ${LOG} \e[0m"

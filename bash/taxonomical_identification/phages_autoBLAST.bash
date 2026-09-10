#!/bin/bash
# Script for automatic BLASTn for bacteriophage taxonomy identification
set -euo pipefail

################
# GLOBAL SETUP #
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/phages"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$(dirname "${SCRIPT_DIR}")")"
TOOL_PATH="${TOOL_PATH:-${REPO_DIR}/python/ncbi}"

# Output directories
BLASTN_NR="${WORK_PATH}/BLASTn_nr_phages"
mkdir -p "${BLASTN_NR}"

# Custom Python scripts
if [ -f "${TOOL_PATH}/ncbi_blastn_auto.py" ]; then
	autoBLASTn="${TOOL_PATH}/ncbi_blastn_auto.py"
elif [ -f "${REPO_DIR}/python/ncbi/ncbi_blastn_auto.py" ]; then
	autoBLASTn="${REPO_DIR}/python/ncbi/ncbi_blastn_auto.py"
else
	autoBLASTn="/storage/student9/tools/ncbi_blastn_auto.py"
fi

if [ -f "${TOOL_PATH}/ncbi_blastn_combine_tsv.py" ]; then
	combineBLASTtsv="${TOOL_PATH}/ncbi_blastn_combine_tsv.py"
elif [ -f "${REPO_DIR}/python/ncbi/ncbi_blastn_combine_tsv.py" ]; then
	combineBLASTtsv="${REPO_DIR}/python/ncbi/ncbi_blastn_combine_tsv.py"
else
	combineBLASTtsv="/storage/student9/tools/ncbi_blastn_combine_tsv.py"
fi

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${WORK_PATH}/phages_autoBLAST.log"
FAIL_LOG="${WORK_PATH}/phages_autoBLAST.failed_skipped.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================" >> "${FAIL_LOG}"
echo " Failure & Skip Log — Started: $(date)" >> "${FAIL_LOG}"
echo "======================================================" >> "${FAIL_LOG}"

log_failure() {
	local phase="$1"
	local sample="$2"
	local status="$3" # e.g. "SKIPPED_EXISTS", "INPUT_MISSING", "EXECUTION_FAILED", "OUTPUT_MISSING"
	local reason="$4"
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${phase}] [${sample}] [${status}] ${reason}" >> "${FAIL_LOG}"
}

# Phage samples ranging from id A21 to A24:
#   WS2762512A21 (E72_BM)
#   WS2762512A22 (TL_1_3_E72)
#   WS2762512A23 (TL_8_E72)
#   WS2762512A24 (E72_TL12_1)

echo "======================================================"
echo " Auto BLASTN for phage samples - Started: $(date)"
echo " Full log           : ${LOG}"
echo " Failure & skip log : ${FAIL_LOG}"
echo "======================================================"

	###########################
	# AUTO BLASTN (NR DATABASE)
	###########################

echo -e "\e[33m ====================================================== \e[0m"
echo -e "\e[33m PHASE: AUTOMATIC NCBI BLASTN SEARCH                   \e[0m"
echo -e "\e[33m ====================================================== \e[0m"

for i in {21..24}
do
	sample_id="WS2762512A${i}"
	BLASTN_NR_SAMPLE="${BLASTN_NR}/${sample_id}"
	expected_tsv="${BLASTN_NR_SAMPLE}/${sample_id}.blastn.combined.tsv"

	echo -e "\e[31m ================================= \e[0m"
	echo -e "\e[31m BLASTN: ${sample_id} (NR)         \e[0m"
	echo -e "\e[31m ================================= \e[0m"

	filtered_assembly="${WORK_PATH}/assembly_phages/${sample_id}.contigs.filtered.fasta"

	if [ ! -f "${filtered_assembly}" ]; then
		echo -e "\e[31m   [${sample_id}] Input assembly missing: ${filtered_assembly} - skipping \e[0m"
		log_failure "Phase_AutoBLASTN" "${sample_id}" "INPUT_MISSING" "File ${filtered_assembly} not found"
		continue
	fi

	if [ -s "${expected_tsv}" ]; then
		echo -e "\e[32m   [${sample_id}] BLASTN results already exist - skipping \e[0m"
		log_failure "Phase_AutoBLASTN" "${sample_id}" "SKIPPED_EXISTS" "Output ${expected_tsv} already exists"
		continue
	fi

	mkdir -p "${BLASTN_NR_SAMPLE}"

	# BLASTN for each node
	if ! conda run -n ncbi python3 "${autoBLASTn}" \
		--delay 15 \
		--input "${filtered_assembly}" \
		--output "${BLASTN_NR_SAMPLE}" \
		--email "huonghm.m23bio@usth.edu.vn"; then
		echo -e "\e[31m   [${sample_id}] ERROR: autoBLASTn failed \e[0m"
		log_failure "Phase_AutoBLASTN" "${sample_id}" "EXECUTION_FAILED" "autoBLASTn returned non-zero"
		continue
	fi

	# Combine all nodes into one .tsv for each sample
	if ! conda run -n ncbi python3 "${combineBLASTtsv}" \
		--input  "${BLASTN_NR_SAMPLE}" \
		--output "${expected_tsv}"; then
		echo -e "\e[31m   [${sample_id}] ERROR: combineBLASTtsv failed \e[0m"
		log_failure "Phase_AutoBLASTN" "${sample_id}" "EXECUTION_FAILED" "combineBLASTtsv returned non-zero"
		continue
	fi

	if [ ! -s "${expected_tsv}" ]; then
		echo -e "\e[31m   [${sample_id}] ERROR: expected combined TSV missing or empty \e[0m"
		log_failure "Phase_AutoBLASTN" "${sample_id}" "OUTPUT_MISSING" "File ${expected_tsv} empty"
		continue
	fi

	echo -e "\e[32m Finish autoBLASTN (nr database) for sample: ${sample_id} \e[0m"
done

echo ""
echo "========================================================"
echo " PHAGE TAXONOMIC IDENTIFICATION COMPLETE - $(date)"
echo "========================================================"
echo "  AutoBLASTN outputs     : ${BLASTN_NR}"
echo ""
echo -e "\e[32m -- LOG FILES ------------------------------------------ \e[0m"
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
echo "========================================================"

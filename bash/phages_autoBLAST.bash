#!/bin/bash
#script for automatical BLASTn for bacteriophage taxonomy identification
#Building No.1
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/phages/"
TOOL_PATH="/storage/student9/tools"

#create working directories
BLASTN_NR="${WORK_PATH}/BLASTn_nr_phages"

mkdir -p "${WORK_PATH}/BLASTn_nr_phages"

#custom .py scripts
autoBLASTn="${TOOL_PATH}/ncbi_blastn_auto.py"
combineBLASTtsv="${TOOL_PATH}/ncbi_blastn_combine_tsv.py"

#create global log
LOG="${WORK_PATH}/phages_autoBLAST.log"
exec > >(tee -a "${LOG}") 2>&1

#phages samples are ranging from id A21 to A24
	#WS2762512A21 (E72_BM)
	#WS2762512A22 (TL_1_3_E72)
	#WS2762512A23 (TL_8_E72)
	#WS2762512A24 (E72_TL12_1)

echo "======================================================"
echo " Auto BLASTN for phage samples - Started: $(date)"
echo "======================================================"

	###########################
	#AUTO BLASTN (NR DATABASE)#
	###########################

	rm -rf "${WORK_PATH}/BLASTn_nr_phages"

for i in {21..24}
do
	sample_id="${i}"

	mkdir -p "${BLASTN_NR}/WS2762512A${sample_id}"
	BLASTN_NR_SAMPLE="${BLASTN_NR}/WS2762512A${sample_id}"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m BLASTN: WS2762512A${sample_id} (NR) \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	filtered_assembly="${WORK_PATH}/assembly_phages/WS2762512A${sample_id}.contigs.filtered.fasta"

	#BLASTN for each node
	python3 ${autoBLASTn} \
	--delay 15 \
	--input "${filtered_assembly}" \
	--output "${BLASTN_NR_SAMPLE}" \
	--email "huonghm.m23bio@usth.edu.vn"

	#combine all nodes into one .tsv for each sample
	python3 ${combineBLASTtsv} \
    	--input  "${BLASTN_NR_SAMPLE}" \
    	--output "${BLASTN_NR_SAMPLE}/WS2762512A${sample_id}.blastn.combined.tsv"

	echo -e "\e[32m Finish autoBLASTN (nr database) for sample: WS2762512A${sample_id} \e[0m"
done

echo "========================================================"
echo " COMPARATIVE GENOMICS COMPLETE - $(date)"
echo "========================================================"
echo "  AutoBLASTN      : ${WORK_PATH}/BLASTn_nr_phages"
echo "========================================================"

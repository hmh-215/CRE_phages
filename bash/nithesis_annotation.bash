#!/bin/bash
#script for annotating bacteriophage ILLUMINA samples
#Building No.1
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/project/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/nithesis"

#Paths to databases
checkvdb="${REF_PATH}/checkv-db-v1.5"
eggnog_db="${REF_PATH}/eggnog_db"

#Create working directories
mkdir -p "${WORK_PATH}/structure_nithesis"
mkdir -p "${WORK_PATH}/annotation_nithesis"
mkdir -p "${ANNOTATION_NITHESIS}/WS2762512A22.term"

STRUCTURE_NITHESIS="${WORK_PATH}/structure_nithesis"
PHAGETERM_REF="${STRUCTURE_NITHESIS}/phageterm_references"
ANNOTATION_NITHESIS="${WORK_PATH}/annotation_nithesis"
PHAGETERM_REF="${WORK_PATH}/phageterm_references"
TERM_NITHESIS="${ANNOTATION_NITHESIS}/WS2762512A22.term"

#Create global log
LOG="${WORK_PATH}/nithesis_annotation.log"
exec > >(tee -a "${LOG}") 2>&1

#annotating nithesis-like phage isolated from samples A22, A23, A24

echo "============================================================"
echo " Nithesis-like phage annotations - Started: $(date)"
echo "============================================================"

	###########################
	#CHECK GENOME COMPLETENESS#
	###########################

for i in {22..24}
do
	sample_id="${i}"

	rm -rf "${STRUCTURE_NITHESIS}/WS2762512A${sample_id}.checkv"
	
	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m CHECKV: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"
	
	filtered_assembly="${SAMPLE_PATH}/phages/assembly_phages/WS2762512A${sample_id}.contigs.filtered.fasta"
	
	#check assembled genome completeness
	conda run -n BPstructure checkv end_to_end \
	-t 16 \
	-d ${checkvdb} \
	"${filtered_assembly}" \
	"${STRUCTURE_NITHESIS}/WS2762512A${sample_id}.checkv"
done

	#############################################
	#PREDICTION OF PHAGE TERMINI & REORIENTATION#
	#############################################

for i in {23..24}
do
	sample_id="${i}"

	mkdir -p "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.term"
	TERM_NITHESIS="${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.term"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m PHAGETERM: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================== \e[0m"
	
	phageterm_ref="${PHAGETERM_REF}/A${sample_id}_NODE_1_metaSPAdes.fasta"
	
	read1f="${WORK_PATH}/mapping_phages/WS2762512A${sample_id}/WS2762512A${sample_id}.R1.f12.fastq.gz"
	read2f="${WORK_PATH}/mapping_phages/WS2762512A${sample_id}/WS2762512A${sample_id}.R2.f12.fastq.gz"

	#find termini of assembled genomes, with reference of the longest contigs (A21, A23, A24)
	conda run -n phageterm_env phageterm \
    	-f "${read1f}" \
    	-p "${read2f}" \
	-r "${phageterm_ref}" \
    	-c 16 \
	-o ${TERM_NITHESIS} \
    	--report_title "WS2762512A${sample_id}_NODE_1.term"

	echo -e "\e[31m Phageterm analysis finished for sample WS2762512A${sample_id}_NODE_1 \e[0m"
done


	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m PHAGETERM: WS2762512A22 \e[0m"
	echo -e "\e[31m ======================= \e[0m"
	
	phageterm_ref="${PHAGETERM_REF}/A22_NODE_2_metaSPAdes.fasta"
	
	read1f="${WORK_PATH}/mapping_phages/WS2762512A22/WS2762512A22.R1.f12.fastq.gz"
	read2f="${WORK_PATH}/mapping_phages/WS2762512A22/WS2762512A22.R2.f12.fastq.gz"

	#find termini of assembled genomes, A22 has 2 different references
	conda run -n phageterm_env ${phageterm} \
    	-f "${read1f}" \
    	-p "${read2f}" \
	-r "${phageterm_ref}" \
    	-c 16 \
	-o ${TERM_NITHESIS} \
    	> "${TERM_NITHESIS}/WS2762512A22.term.out" \
	2> "${TERM_NITHESIS}/WS2762512A22.term.log"

	echo "done phageterm"

#Extract each nodes of assemblies of above 90% completeness
for i in {21..24}
do
	sample_id="${i}"
	sample_name="WS2762512A${sample_id}"

	checkv_tsv="${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.checkv/quality_summary.tsv"
	filtered_assembly="${SAMPLE_PATH}/phages/assembly_phages/WS2762512A${sample_id}.contigs.filtered.fasta"

	if [[ ! -f "${checkv_tsv}" ]]; then
		echo "WARNING: CheckV summary not found for ${sample_name}, skipping node extraction."
		continue
	fi

	# Parse TSV: column 1 = contig_id, column 9 = completeness (%)
	# Skip header line; filter rows where completeness (col 9) >= 90
	tail -n +2 "${checkv_tsv}" | awk -F'\t' '$9 != "NA" && $9+0 >= 90 { print $1 }' | \

	while read -r contig_id
	do
		# Extract the NODE number from contig_id (e.g. "NODE_1_length_..." -> "1")
		node_num=$(echo "${contig_id}" | grep -oP '(?<=NODE_)\d+')

		if [[ -z "${node_num}" ]]; then
			echo "WARNING: Could not parse node number from contig ID '${contig_id}' in ${sample_name}, skipping."
			continue
		fi

		out_fasta="${PHAGETERM_REF}/${sample_name}_NODE_${node_num}.metaSPAdes.fasta"
		new_header="${sample_name}_NODE_${node_num}.metaSPAdes.fasta"

		# Extract the single contig and rename the header
		conda run -n assembly seqkit grep \
		-p "${contig_id}" \
		"${filtered_assembly}" \
		| conda run -n assembly seqkit replace \
		--pattern ".*" \
		--replacement "${new_header}" \
		> "${out_fasta}"

		echo "  Extracted: ${out_fasta}  (contig: ${contig_id})"
	done
done

	echo -e "\e[32m Node extraction complete. Files written to: ${PHAGETERM_REF} \e[0m"

for phageterm_ref in "${PHAGETERM_REF}"/*.metaSPAdes.fasta
do
	[[ -f "${phageterm_ref}" ]] || continue

	# Parse sample name and node from filename, e.g. WS2762512A21_NODE_1.metaSPAdes.fasta
	basename_no_ext=$(basename "${phageterm_ref}" .metaSPAdes.fasta)  # WS2762512A21_NODE_1
	# sample_name = everything before the last _NODE_N suffix
	sample_name=$(echo "${basename_no_ext}" | grep -oP '^WS\d+A\d+')
	node_num=$(echo "${basename_no_ext}" | grep -oP '(?<=NODE_)\d+')

	if [[ -z "${sample_name}" || -z "${node_num}" ]]; then
		echo "WARNING: Could not parse sample/node from '${phageterm_ref}', skipping."
		continue
	fi

	TERM_NITHESIS="${STRUCTURE_NITHESISS}/${sample_name}_NODE_${node_num}.term"
	mkdir -p "${TERM_NITHESIS}"
	cd "@{TERM_NITHESIS}"

	read1f="${MAPPING_PHAGES}/WS2762512A${sample_id}/WS2762512A${sample_id}.R1.f12.fastq.gz"
	read2f="${MAPPING_PHAGES}/WS2762512A${sample_id}/WS2762512A${sample_id}.R2.f12.fastq.gz"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m PHAGETERM: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	#find termini of assembled genomes
	conda run -n phageterm_env phageterm \
    	-f "${read1f}" \
    	-p "${read2f}" \
	-r "${phageterm_ref}" \
    	-c ${threads} \
	--report_title "WS2762512A${sample_id}_NODE_${node_num}.term"

	echo -e "\e[32m PhageTerm analysis finish for sample WS2762512A${sample_id} \e[0m"
done


	#############################
	#MULTIPLE SEQUENCE ALIGNMENT#
	#############################

	mafft_input="${ANNOTATION_NITHESIS}/mafft_input_nithesis.fasta"

	conda run -n BPannotate mafft \
	--thread 16 \
	--clustalout \
	"${mafft_input}" > "${ANNOTATION_NITHESIS}/mafft_output_nithesis.aln"

	#############
	#ANNOTATIONS#
	#############

#annotations for sample A22

	sample_id="22"	

	mkdir -p "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.pharokka"

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m PHAROKKA: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	assembly_nithesis="${PHAGETERM_REF}/A22_NODE_2_metaSPAdes.fasta"

	#annotation with prodigal (gold standard method)
	conda run -n BPannotate pharokka.py \
	-t 16 -f \
	-d "${REF_PATH}/pharokka_db" \
	-i "${assembly_nithesis}" \
	-g prodigal \
	-o "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.pharokka" \
	-p "WS2762512A${sample_id}.pharokka"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m PHANOTATE: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#additional annotation for comparision with pharokka (alternative ORF caller)
	conda run -n BPannotate phanotate.py \
	--format tabular \
	-o "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.phanotate" \
	"${assembly_nithesis}"

	mkdir -p "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.eggnog"
	EGGNOG_NITHESIS="${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.eggnog"

#annotate for A23, A24

for i in {23..24}
do
	sample_id="${i}"
	
	mkdir -p "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.pharokka"

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m PHAROKKA: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	assembly_nithesis="${PHAGETERM_REF}/A${sample_id}_NODE_1_metaSPAdes.fasta"

	#annotation with prodigal (gold standard method)
	conda run -n BPannotate pharokka.py \
	-t 16 -f \
	-d "${REF_PATH}/pharokka_db" \
	-i "${assembly_nithesis}" \
	-g prodigal \
	-o "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.pharokka" \
	-p "WS2762512A${sample_id}.pharokka"

done

for i in {23..24}
do
	sample_id="${i}"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m PHANOTATE: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	assembly_nithesis="${PHAGETERM_REF}/A${sample_id}_NODE_1_metaSPAdes.fasta"

	#additional annotation for comparision with pharokka (alternative ORF caller)
	conda run -n BPannotate phanotate.py \
	--format tabular \
	-o "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.phanotate" \
	"${assembly_nithesis}"
done

for i in {22..24}
do
	sample_id="${i}"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EGGNOG-MAPPER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	mkdir -p "${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.eggnog"
	EGGNOG_NITHESIS="${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.eggnog"

	prodigal_input="${ANNOTATION_NITHESIS}/WS2762512A${sample_id}.pharokka/prodigal.faa"

	#annotate GO terms
	# Step 1: DIAMOND search only
	conda run -n BPannotate emapper.py \
	--cpu 16 \
	--dbmem \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_NITHESIS}" \
	--output_dir "${EGGNOG_NITHESIS}" \
	--no_annot \
	--override \
	-i "${prodigal_input}" \
	-o "WS2762512A${sample_id}.eggnog"

	# Step 2: Annotation only using hits from step 1
	conda run -n BPannotate emapper.py \
	--cpu 16 \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_NITHESIS}" \
	--output_dir "${EGGNOG_NITHESIS}" \
	--annotate_hits_table "${EGGNOG_NITHESIS}/WS2762512A${sample_id}.eggnog.emapper.seed_orthologs" \
	-m no_search \
	--excel \
	--tax_scope Viruses \
	--override \
	-o "WS2762512A${sample_id}.eggnog"
	
	echo -e "\e[32m EggNOG-mapper analysis finish for sample WS2762512A${sample_id} \e[0m"
done

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - $(date)"
echo "========================================================"
echo "  Annotation          : ${ANNOTATION_NITHESIS}"
echo "  Genome completness  : ${STRUCTURE_NITHESIS}/*.checkv"
echo "  Termini             : ${TERM_NITHESIS}"
echo "  Full log            : ${LOG}"
echo "========================================================"

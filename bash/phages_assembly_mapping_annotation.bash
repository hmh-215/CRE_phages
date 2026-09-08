#!/bin/bash
#script for mapping, assembly, and annotations of bacteriophage ILLUMINA samples
#Building No.1
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/phages"
threads="16"

#path to databases and references
E72_reference="${SAMPLE_PATH}/Escherichia_coli/assembly_Ecoli/WS2762512A08.contigs.filtered.fasta"

checkvdb="${REF_PATH}/checkv-db-v1.5"
eggnog_db="${REF_PATH}/eggnog_db"

truseq2="${REF_PATH}/TruSeq2-PE.fa"
truseq3="${REF_PATH}/TruSeq3-PE.fa"

#Create working directories
PREPROCESSING_PHAGES="${WORK_PATH}/preprocessing_phages"
MAPPING_PHAGES="${WORK_PATH}/mapping_phages"
ASSEMBLY_PHAGES="${WORK_PATH}/assembly_phages"
QUAST_RAW_PHAGES_INPUTS="${ASSEMBLY_PHAGES}/quast_raw_phages_inputs"
STRUCTURE_PHAGES="${WORK_PATH}/structure_phages"
PHAGETERM="${STRUCTURE_PHAGES}/phageterm"
PHAGETERM_REF="${PHAGETERM}/phageterm_references"
ANNOTATION_PHAGES="${WORK_PATH}/annotation_phages"

mkdir -p "${WORK_PATH}/preprocessing_phages"
mkdir -p "${WORK_PATH}/mapping_phages"
mkdir -p "${PREPROCESSING_PHAGES}/fastqc_phages_raw"
mkdir -p "${PREPROCESSING_PHAGES}/fastqc_phages_pp"
mkdir -p "${MAPPING_PHAGES}/fastqc_phages_f12"
mkdir -p "${WORK_PATH}/assembly_phages"
mkdir -p "${ASSEMBLY_PHAGES}/quast_raw_phages_inputs"
mkdir -p "${WORK_PATH}/structure_phages"
mkdir -p "${STRUCTURE_PHAGES}/phageterm"
mkdir -p "${PHAGETERM}/phageterm_references"
mkdir -p "${WORK_PATH}/annotation_phages"

#create global log
LOG="${WORK_PATH}/phages_assembly_mapping_annotation.log"
exec > >(tee -a "${LOG}") 2>&1

#phages samples are ranging from id A21 to A24
	#WS2762512A21 (E72_BM)
	#WS2762512A22 (TL_1_3_E72)
	#WS2762512A23 (TL_8_E72)
	#WS2762512A24 (E72_TL12_1) 
#However, there are bacterial DNA involved that must be removed.

echo "============================================================"
echo " Phages assembly, mapping, and annotations - Started: $(date)"
echo "============================================================"

	#########################
	#SAMPLE DATA PREPARATION#
	#########################

for i in {21..24}
do
	sample_id="${i}"

	read1="${SAMPLE_PATH}/WS2762512A${sample_id}_R1.fastq.gz"
	read2="${SAMPLE_PATH}/WS2762512A${sample_id}_R2.fastq.gz"
	trim_output_1="${PREPROCESSING_PHAGES}/WS2762512A${sample_id}.R1"
	trim_output_2="${PREPROCESSING_PHAGES}/WS2762512A${sample_id}.R2"
	read1t="${PREPROCESSING_PHAGES}/WS2762512A${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_PHAGES}/WS2762512A${sample_id}.R2.paired.fastq.gz"

	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: WS2762512A${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	#review raw fastq
	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${PREPROCESSING_PHAGES}/fastqc_phages_raw" \
	"${read1}" "${read2}"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m TRIMMOMATIC: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#samples data manipulation
	conda run -n preprocessing trimmomatic PE \
	-threads ${threads} \
	-phred64 \
	 ${read1} ${read2} \
	${trim_output_1}.paired.fastq.gz ${trim_output_1}.unpaired.fastq.gz \
	${trim_output_2}.paired.fastq.gz ${trim_output_2}.unpaired.fastq.gz \
	ILLUMINACLIP:"${truseq2}":2:30:10:2:True \
	ILLUMINACLIP:"${truseq3}":2:30:10:2:True \
	LEADING:3 TRAILING:3 MINLEN:50

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m FASTQC: WS2762512A${sample_id} (PP) \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#review trimmed fastq
	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${PREPROCESSING_PHAGES}/fastqc_phages_pp" \
	"${read1t}" "${read2t}"
done

	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m MULTIQC: PHAGE SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	#review all raw fastq
	conda run -n preprocessing multiqc \
	--force \
	"${PREPROCESSING_PHAGES}/fastqc_phages_raw" \
	--filename "${PREPROCESSING_PHAGES}/multiqc_phages_raw"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m MULTIQC: PHAGE SAMPLES (PP) \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	#review all trimmed fastq
	conda run -n preprocessing multiqc \
	--force \
	"${PREPROCESSING_PHAGES}/fastqc_phages_pp" \
	--filename "${PREPROCESSING_PHAGES}/multiqc_phages_pp"

	####################################
	#MAPPING & FILTERING UNMAPPED READS#
	####################################

	conda run -n mapping bwa index "${E72_reference}"

for i in {21..24}
do
	sample_id="${i}"

	mkdir -p "${MAPPING_PHAGES}/WS2762512A${sample_id}"

	MAPPING_PHAGE_SAMPLE="${MAPPING_PHAGES}/WS2762512A${sample_id}"
	read1t="${PREPROCESSING_PHAGES}/WS2762512A${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_PHAGES}/WS2762512A${sample_id}.R2.paired.fastq.gz"
	mapped_bam="${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.bam"
	sorted_bam="${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sorted_coord.bam"
	read1f="${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R1.f12.fastq.gz"
	read2f="${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R2.f12.fastq.gz"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m BWA MEM: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#Map trimmed fastq to E72 references
	conda run -n mapping bwa mem \
	-t ${threads} \
	"${E72_reference}" \
	"${read1t}" "${read2t}" \
	> "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sam"

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m SAMTOOLS: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	#Conversion from .sam to .bam format
	conda run -n mapping samtools view \
	-@ ${threads} \
	-bS \
	-o "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.bam" \
	"${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sam"

	rm "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sam"

	#sorting mapped reads based on coordinates
	conda run -n mapping samtools sort \
	-@ ${threads} \
	-o "${sorted_bam}" \
	"${mapped_bam}"

	#indexing mapped reads
	conda run -n mapping samtools index \
	-@ ${threads} \
	"${sorted_bam}"

	#examine mapped reads coverage
	conda run -n mapping samtools coverage \
	"${sorted_bam}" \
	> "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sorted_coord.coverage.txt"

	#filter unmapped reads
	conda run -n mapping samtools view \
	-@ ${threads} \
	-f 12 \
	-b \
	"${sorted_bam}" \
	-o  "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.f12.bam"
	#for online conversion of bitwise flags, see: https://42basepairs.com/tools/sam-flag

	#sort queries by names
	conda run -n mapping samtools sort \
	-@ ${threads} \
	-n \
	-o "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sorted_name.f12.bam" \
	"${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.f12.bam"

	#converting bam to fastq
	conda run -n mapping samtools bam2fq \
	-@ ${threads} -n \
	-1 "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R1.f12.fastq.gz" \
	-2 "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R2.f12.fastq.gz" \
	"${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.sorted_name.f12.bam"

	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m FASTQC: WS2762512A${sample_id} (-f 12) \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	#review trimmed fastq
	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${MAPPING_PHAGES}/fastqc_phages_f12" \
	"${read1f}" "${read2f}"
done

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m MULTIQC: PHAGE SAMPLES (-f 12) \e[0m"
	echo -e "\e[31m ============================== \e[0m"

	#review all trimmed fastq
	conda run -n preprocessing multiqc \
	--force \
	"${MAPPING_PHAGES}/fastqc_phages_f12" \
	--filename "${MAPPING_PHAGES}/multiqc_phages_f12"

	##################
	#DE NOVO ASSEMBLY#
	##################

for i in {21..24}
do
	sample_id="${i}"

	MAPPING_PHAGE_SAMPLE="${MAPPING_PHAGES}/WS2762512A${sample_id}"
	read1f="${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R1.f12.fastq.gz"
	read2f="${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R2.f12.fastq.gz"
	assembly="${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.contig.fasta/contigs.fasta"
	filtered_assembly="${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.contigs.filtered.fasta"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SEQKIT: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#checking coverage of fastq files
	conda run -n assembly seqkit stats \
	--threads ${threads} \
	--all --tabular \
	-o "${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.f12.fastq.stats.txt"\
	"${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R1.f12.fastq.gz" \
	"${MAPPING_PHAGE_SAMPLE}/WS2762512A${sample_id}.R2.f12.fastq.gz"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SPADES: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#de novo assembly phage contigs
	conda run -n assembly spades.py \
	--meta \
	-1 "${read1f}" \
	-2 "${read2f}" \
	-t ${threads} \
	-o "${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.contig.fasta"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#review raw assembly statistics
	conda run -n assembly quast.py \
	"${assembly}" \
	--threads ${threads} \
	--output-dir "${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.raw.quast"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SEQKIT: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#remove contigs with less than 500bp
	conda run -n assembly seqkit seq \
	--min-len 500 \
	"${assembly}" \
	> "${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.contigs.filtered.fasta"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} (PP) \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	#review filtered assembly statistics
	conda run -n assembly quast.py \
	"${filtered_assembly}" \
	--threads ${threads} \
	--output-dir "${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.filtered.quast"
done

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m QUAST: ALL PHAGE SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ============================== \e[0m"

for i in {21..24}
do
	sample_id="${i}"

	cp -f "${ASSEMBLY_PHAGES}/WS2762512A${sample_id}.contig.fasta/contigs.fasta" \
	"${QUAST_RAW_PHAGES_INPUTS}/WS2762512A${sample_id}.contigs.raw.fasta"
done

	#review all raw assemblies statistics
	conda run -n assembly quast.py \
	--threads ${threads} \
	--output-dir "${ASSEMBLY_PHAGES}/all_phages.raw.quast" \
	"${QUAST_RAW_PHAGES_INPUTS}/"*.contigs.raw.fasta

	echo -e "\e[31m ============================= \e[0m"
	echo -e "\e[31m QUAST: ALL PHAGE SAMPLES (PP) \e[0m"
	echo -e "\e[31m ============================= \e[0m"

	#review all filtered assemblies statistics
	conda run -n assembly quast.py \
	--threads ${threads} \
	--output-dir "${ASSEMBLY_PHAGES}/all_phages.raw.quast" \
	"${ASSEMBLY_PHAGES}/"*.contigs.filtered.fasta

	###########################
	#CHECK GENOME COMPLETENESS#
	###########################

for i in {21..24}
do
	sample_id="${i}"

	rm -rf "${STRUCTURE_PHAGES}/WS2762512A${sample_id}.checkv"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m CHECKV: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	filtered_assembly="${WORK_PATH}/assembly_phages/WS2762512A${sample_id}.contigs.filtered.fasta"

	#check assembled genome completeness
	conda run -n BPstructure checkv end_to_end \
	-t ${threads} \
	-d ${checkvdb} \
	"${filtered_assembly}" \
	"${STRUCTURE_PHAGES}/WS2762512A${sample_id}.checkv"
done

	#############################################
	#PREDICTION OF PHAGE TERMINI & REORIENTATION#
	#############################################

	#Step 1 — Extract contigs with CheckV completeness >= 90%
	echo -e "\e[32m Extracting high-completeness contigs (>=90%) from CheckV results... \e[0m"

for i in {21..24}
do
	sample_id="${i}"
	sample_name="WS2762512A${sample_id}"

	checkv_tsv="${STRUCTURE_PHAGES}/${sample_name}.checkv/quality_summary.tsv"
	filtered_assembly="${ASSEMBLY_PHAGES}/${sample_name}.contigs.filtered.fasta"

	#Guard: skip if CheckV output missing
	if [[ ! -f "${checkv_tsv}" ]]; then
		echo "WARNING: CheckV summary not found for ${sample_name}, skipping."
		continue
	fi

	#Guard: skip if assembly missing
	if [[ ! -f "${filtered_assembly}" ]]; then
		echo "WARNING: Filtered assembly not found for ${sample_name}, skipping."
		continue
	fi

	echo -e "\e[32m Processing CheckV results for ${sample_name}... \e[0m"

	#Parse CheckV TSV:
	#Column 1 = contig_id
	#Column 9 = completeness (%) — skip NA values, keep >= 90%
	#tail -n +2 skips the header line
	tail -n +2 "${checkv_tsv}" \
	| awk -F'\t' '$10 != "NA" && $10+0 >= 90 { print $1 }' \
	| while read -r contig_id
	do
		#Extract NODE number from SPAdes contig ID
		#e.g. NODE_1_length_72827_cov_610.9 -> 1
		node_num=$(echo "${contig_id}" | grep -oP '(?<=NODE_)\d+')

		if [[ -z "${node_num}" ]]; then
			echo "  WARNING: Could not parse NODE number from '${contig_id}', skipping."
			continue
		fi

		out_fasta="${PHAGETERM_REF}/${sample_name}_NODE_${node_num}.metaSPAdes.fasta"

		#Use a clean header — just sample + node, no .fasta extension in header
		new_header="${sample_name}_NODE_${node_num}"

		echo "  Extracting contig: ${contig_id} ? ${out_fasta}"

		#Extract the matching contig and rename its header
		tmp_fa="${PHAGETERM_REF}/tmp.${sample_name}.${node_num}.fa"

		conda run -n assembly seqkit grep \
		-p "${contig_id}" \
		"${filtered_assembly}" \
		> "${tmp_fa}"

		conda run -n assembly seqkit replace \
		-p '^.*$' \
		-r "${new_header}" \
		"${tmp_fa}" \
		> "${out_fasta}"

		rm -f "${tmp_fa}"


		#Verify extraction produced a non-empty file
		if [[ ! -s "${out_fasta}" ]]; then
			echo "  WARNING: Extracted file is empty for ${contig_id} — check contig ID matches assembly headers."
			rm -f "${out_fasta}"
			continue
		fi

		echo -e "\e[32m Saved: ${out_fasta} \e[0m"
	done

	#Report how many contigs were extracted for this sample
	extracted=$(ls "${PHAGETERM_REF}/${sample_name}"_NODE_*.metaSPAdes.fasta 2>/dev/null | wc -l)
	echo -e "\e[32m ${sample_name}: ${extracted} contig(s) extracted with >=90% completeness \e[0m"
done

	echo -e "\e[32m Node extraction complete. Files written to: ${PHAGETERM_REF} \e[0m"
	ls -lh "${PHAGETERM_REF}/"

	echo -e "\e[32m Starting PhageTerm analysis on extracted contigs... \e[0m"

	#change working directory to phageterm
	cd "${PHAGETERM}"

for phageterm_ref in "${PHAGETERM_REF}"/*.metaSPAdes.fasta
do
	#Skip if glob matched nothing
	[[ -f "${phageterm_ref}" ]] || continue

	#Parse sample_name and node_num from filename
	#e.g. WS2762512A21_NODE_1.metaSPAdes.fasta
	basename_no_ext=$(basename "${phageterm_ref}" .metaSPAdes.fasta)
	sample_name=$(echo "${basename_no_ext}" | grep -oP '^WS\d+A\d+')
	sample_id=${sample_name#WS2762512A}
	node_num=$(echo "${basename_no_ext}"    | grep -oP '(?<=NODE_)\d+')

	#Guard: skip if parsing failed
	if [[ -z "${sample_name}" || -z "${node_num}" ]]; then
		echo -e "\e[32m WARNING: Could not parse sample/node from '${phageterm_ref}', skipping. \e[0m"
		continue
	fi

	#BUG FIX 2: derive read paths from sample_name, NOT from stale ${sample_id}
	#sample_name is correctly parsed from the filename above
	read1f="${MAPPING_PHAGES}/${sample_name}/${sample_name}.R1.f12.fastq.gz"
	read2f="${MAPPING_PHAGES}/${sample_name}/${sample_name}.R2.f12.fastq.gz"

	#Guard: skip if reads are missing
	if [[ ! -f "${read1f}" || ! -f "${read2f}" ]]; then
		echo -e "\e[32m WARNING: Filtered reads not found for ${sample_name}: \e[0m"
		echo -e "\e[32m R1: ${read1f} \e[0m"
		echo -e "\e[32m R2: ${read2f} \e[0m"
		echo -e "\e[32m Skipping PhageTerm for this contig. \e[0m"
		continue
	fi

	#BUG FIX 1 + 3: create dedicated output directory and pass via --DR_path
	TERM_PHAGE="${STRUCTURE_PHAGES}/${sample_name}_NODE_${node_num}.term"
	mkdir -p "${TERM_PHAGE}"

	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m PHAGETERM: ${sample_name} NODE_${node_num} \e[0m"
	echo -e "\e[31m ================================ \e[0m"

	echo -e "\e[32m Reference : ${phageterm_ref} \e[0m"
	echo -e "\e[32m Reads R1  : ${read1f} \e[0m"
	echo -e "\e[32m Reads R2  : ${read2f} \e[0m"
	echo -e "\e[32m Output dir: ${TERM_PHAGE} \e[0m"

	conda run -n phageterm_env phageterm \
	-f "${read1f}" \
	-p "${read2f}" \
	-r "${phageterm_ref}" \
	-c "${threads}" \
	--report_title "${sample_name}_NODE_${node_num}"

	echo -e "\e[32m PhageTerm complete: ${sample_name} NODE_${node_num} \e[0m"
	echo -e "\e[32m Output files in: ${TERM_PHAGE} \e[0m"
	ls -lh "${TERM_PHAGE}/"
done

	echo -e "\e[32m All PhageTerm analyses complete. \e[0m"
	echo -e "\e[32m Results stored under: ${STRUCTURE_PHAGES}/ \e[0m"

	#############
	#ANNOTATIONS#
	#############

for i in {21..24}
do
	sample_id="${i}"

	mkdir -p "${ANNOTATION_PHAGES}/WS2762512A${sample_id}.pharokka"
	mkdir -p "${ANNOTATION_PHAGES}/WS2762512A${sample_id}.eggnog"

	EGGNOG_PHAGE="${ANNOTATION_PHAGES}/WS2762512A${sample_id}.eggnog"
	filtered_assembly="${WORK_PATH}/assembly_phages/WS2762512A${sample_id}.contigs.filtered.fasta"
	prodigal_input="${ANNOTATION_PHAGES}/WS2762512A${sample_id}.pharokka/prodigal.faa"	

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m PHAROKKA: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	#annotation with prodigal (gold standard method)
	conda run -n BPannotation pharokka.py \
	-t ${threads} -f \
	-d "${REF_PATH}/pharokka_db" \
	-i "${filtered_assembly}" \
	-g prodigal \
	-o "${ANNOTATION_PHAGES}/WS2762512A${sample_id}.pharokka" \
	-p "WS2762512A${sample_id}.pharokka"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m PHANOTATE: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#additional annotation for comparision with pharokka (alternative ORF caller)
	conda run -n BPannotation phanotate.py \
	--format tabular \
	-o "${ANNOTATION_PHAGES}/WS2762512A${sample_id}.phanotate" \
	${filtered_assembly}

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EGGNOG-MAPPER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	#annotate GO terms
	#Step 1: DIAMOND search only
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--dbmem \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_PHAGE}" \
	--output_dir "${EGGNOG_PHAGE}" \
	--no_annot \
	--override \
	-i "${prodigal_input}" \
	-o "WS2762512A${sample_id}.eggnog"

	#Step 2: Annotation only using hits from step 1
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_PHAGE}" \
	--output_dir "${EGGNOG_PHAGE}" \
	--annotate_hits_table "${EGGNOG_PHAGE}/WS2762512A${sample_id}.eggnog.emapper.seed_orthologs" \
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
echo "  Preprocessing       : ${PREPROCESSING_PHAGES}"
echo "  Mapping             : ${MAPPING_PHAGES}"
echo "  Phage assemblies    : ${ASSEMBLY_PHAGES}"
echo "  Annotation          : ${ANNOTATION_PHAGES}"
echo "  Genome completness  : ${STRUCTURE_PHAGES}/*.checkv"
echo "  Termini             : ${PHAGETERM}"
echo "  Full log            : ${LOG}"
echo "========================================================"
#!/bin/bash
#scripts for de novo assembly and annotation of E. coli ILLUMINA samples
#Building No.1
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/Escherichia_coli"

threads=16

#path to databases
truseq2="${REF_PATH}/TruSeq2-PE.fa"
truseq3="${REF_PATH}/TruSeq3-PE.fa"

bakta_db="${REF_PATH}/bakta_db/db"
amrfinder_db="${REF_PATH}/bakta_db/db/amrfinderplus-db/latest"
eggnog_db="${REF_PATH}/eggnog_db"

#create working directories

PREPROCESSING_ECOLI="${WORK_PATH}/preprocessing_Ecoli"
ASSEMBLY_ECOLI="${WORK_PATH}/assembly_Ecoli"
ANNOTATION_ECOLI="${WORK_PATH}/annotation_Ecoli"
ANNOTATION_ECOLI_CHECKM="${ANNOTATION_ECOLI}/checkm"
QUAST_RAW_ECOLI_INPUTS="${ASSEMBLY_ECOLI}/quast_raw_Ecoli_inputs"
CHECKM_INPUTS="${ANNOTATION_ECOLI}/checkm/inputs_Ecoli"
ANNOTATION_ECOLI_MLST="${ANNOTATION_ECOLI}/mlst"

mkdir -p "${WORK_PATH}"
mkdir -p "${WORK_PATH}/preprocessing_Ecoli"
mkdir -p "${WORK_PATH}/assembly_Ecoli"
mkdir -p "${WORK_PATH}/annotation_Ecoli"
mkdir -p "${ASSEMBLY_ECOLI}/quast_raw_Ecoli_inputs"
mkdir -p "${ANNOTATION_ECOLI}/checkm"
mkdir -p "${ANNOTATION_ECOLI}/checkm/inputs_Ecoli"
mkdir -p "${ANNOTATION_ECOLI}/mlst"
mkdir -p "${ANNOTATION_ECOLI}/bakta"
mkdir -p "${ANNOTATION_ECOLI}/amrfinder"
mkdir -p "${ANNOTATION_ECOLI}/eggnog"

#create global log
LOG="${WORK_PATH}/Ecoli_assembly_annotation.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================"
echo " E. coli assemblies and annotations - Started: $(date)"
echo "======================================================"

	#########################
	#SAMPLE DATA PREPARATION#
	#########################

	mkdir -p "${PREPROCESSING_ECOLI}/fastqc_Ecoli_raw"
	mkdir -p "${PREPROCESSING_ECOLI}/fastqc_Ecoli_pp"

#A08 sample is the host of phages (E. coli)		
for sample_id in "05" "08" "09" "11" "12" "13" "16" "17" "18" "19" "20"
do
	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: WS2762512A${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"
	
	read1="${SAMPLE_PATH}/WS2762512A${sample_id}_R1.fastq.gz"
	read2="${SAMPLE_PATH}/WS2762512A${sample_id}_R2.fastq.gz"

	#review raw fastq
	conda run -n preprocessing fastqc \
   	--threads ${threads} \
    	--outdir "${PREPROCESSING_ECOLI}/fastqc_Ecoli_raw" \
    	"${read1}" "${read2}"	

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m TRIMMOMATIC: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ========================= \e[0m"
	
	trim_output_1="${PREPROCESSING_ECOLI}/WS2762512A${sample_id}.R1"
	trim_output_2="${PREPROCESSING_ECOLI}/WS2762512A${sample_id}.R2"

	#fastq data manipulation	
	conda run -n preprocessing trimmomatic PE \
	-threads ${threads} \
	-phred64 \
	 ${read1} ${read2} \
	${trim_output_1}.paired.fastq.gz ${trim_output_1}.unpaired.fastq.gz \
	${trim_output_2}.paired.fastq.gz ${trim_output_2}.unpaired.fastq.gz \
	ILLUMINACLIP:"${truseq2}":2:30:10:2:True \
	ILLUMINACLIP:"${truseq3}":2:30:10:2:True \
	LEADING:3 TRAILING:3 MINLEN:50
	
	read1t=${trim_output_1}.paired.fastq.gz
	read2t=${trim_output_2}.paired.fastq.gz
	
	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m FASTQC: WS2762512A${sample_id} (PP) \e[0m"
	echo -e "\e[31m ========================= \e[0m"
	
	#review filtered fastq
	conda run -n preprocessing fastqc \
    	--threads ${threads} \
    	--outdir "${PREPROCESSING_ECOLI}/fastqc_Ecoli_pp" \
	"${read1t}" "${read2t}"
done

	echo -e "\e[31m ================================== \e[0m"
	echo -e "\e[31m MULTIQC: ALL E. COLI SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ================================== \e[0m"

	#review all trimmed fastq
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_ECOLI}/fastqc_Ecoli_raw" \
    	--filename "${PREPROCESSING_ECOLI}/multiqc_Ecoli_raw"

	echo -e "\e[31m ================================= \e[0m"
	echo -e "\e[31m MULTIQC: ALL E. COLI SAMPLES (PP) \e[0m"
	echo -e "\e[31m ================================= \e[0m"

	#review all trimmed fastq
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_ECOLI}/fastqc_Ecoli_pp" \
    	--filename "${PREPROCESSING_ECOLI}/multiqc_Ecoli_pp"

	##################
	#DE NOVO ASSEMBLY#
	##################

for sample_id in "05" "08" "09" "11" "12" "13" "16" "17" "18" "19" "20"
do
	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SPADES: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"
	
	read1t="${PREPROCESSING_ECOLI}/WS2762512A${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_ECOLI}/WS2762512A${sample_id}.R2.paired.fastq.gz"

	#de novo assembly
	conda run -n assembly spades.py \
    	-1 "${read1t}" \
    	-2 "${read2t}" \
    	--careful \
    	-t ${threads} \
    	-o "${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.contig.fasta"

	Ecoli_assembly="${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.contig.fasta/contigs.fasta"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================= \e[0m"
	
	#review raw assembly statistics
	conda run -n assembly quast.py \
    	"${Ecoli_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.raw.quast"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SEQKIT: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"
	
	#remove contigs with less than 500bp
	conda run -n assembly seqkit seq \
    	--min-len 500 \
    	"${Ecoli_assembly}" \
    	> "${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.contigs.filtered.fasta"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} (PP) \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	filtered_Ecoli_assembly="${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.contigs.filtered.fasta"

	#review filtered assembly statistics
	conda run -n assembly quast.py \
    	"${filtered_Ecoli_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.filtered.quast"
done

	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m QUAST: ALL E. COLI SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ================================ \e[0m"

for sample_id in "05" "08" "09" "11" "12" "13" "16" "17" "18" "19" "20"
do
	cp -f "${ASSEMBLY_ECOLI}/WS2762512A${sample_id}.contig.fasta/contigs.fasta" "${QUAST_RAW_ECOLI_INPUTS}/WS2762512A${sample_id}.contigs.raw.fasta"
done
	#review all raw assemblies statistics
	conda run -n assembly quast.py \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_ECOLI}/all_Ecoli.raw.quast" \
	"${QUAST_RAW_ECOLI_INPUTS}/"*.contigs.raw.fasta	

	echo -e "\e[31m =============================== \e[0m"
	echo -e "\e[31m QUAST: ALL E. COLI SAMPLES (PP) \e[0m"
	echo -e "\e[31m =============================== \e[0m"

	#review all filtered assemblies statistics
	conda run -n assembly quast.py \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_ECOLI}/all_Ecoli.filtered.quast" \
	"${ASSEMBLY_ECOLI}/"*contigs.filtered.fasta

	###########################
	#CHECK GENOME COMPLETENESS#
	###########################
	
	cp -f "${ASSEMBLY_ECOLI}/"*contigs.filtered.fasta "${CHECKM_INPUTS}"

	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m CHECKM: ALL E. COLI SAMPLES (PP) \e[0m"
	echo -e "\e[31m ================================ \e[0m"

	#check assembly genomes completeness
	conda run -n BPstructure checkm lineage_wf \
	--threads ${threads} \
	--extension fasta \
	--tab_table \
	${CHECKM_INPUTS} \
	${ANNOTATION_ECOLI_CHECKM}

	# CheckM quality assessment and summary
	if [ -f "${ANNOTATION_ECOLI_CHECKM}/lineage.ms" ]; then
		conda run -n BPstructure checkm qa \
		"${ANNOTATION_ECOLI}/checkm/lineage.ms" \
		"${ANNOTATION_ECOLI_CHECKM}" \
		-o 2 \
		--tab_table \
		-f "${ANNOTATION_ECOLI_CHECKM}/all_Ecoli.quality.checkm.tsv"
		echo -e "\e[32m CheckM complete: ${ANNOTATION_ECOLI_CHECKM}/all_Ecoli.quality.checkm.tsv \e[0m"
	else
		echo -e "\e[32m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
		echo -e "\e[32m Genome quality assessed by assembly statistics instead \e[0m"
	fi

	#############################
	#MULTI-LOCUS SEQUENCE TYPING#
	#############################

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m MLST: ALL E. COLI SAMPLES \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#serotyping strains	
	conda run -n BPtyping mlst \
	--threads ${threads} \
	--scheme ecoli \
	--csv \
	"${ASSEMBLY_ECOLI}/"*contigs.filtered.fasta \
	> "${ANNOTATION_ECOLI_MLST}/all_Ecoli.mlst.csv"

	#############
	#ANNOTATIONS#
	#############

for sample_id in "05" "08" "09" "11" "12" "13" "16" "17" "18" "19" "20"
do
	mkdir -p "${ANNOTATION_ECOLI}/bakta/WS2762512A${sample_id}.bakta"
	mkdir -p "${ANNOTATION_ECOLI}/amrfinder/WS2762512A${sample_id}.amrfinder"
	mkdir -p "${ANNOTATION_ECOLI}/eggnog/WS2762512A${sample_id}.eggnog"

	ANNOTATION_ECOLI_BAKTA="${ANNOTATION_ECOLI}/bakta/WS2762512A${sample_id}.bakta"
	ANNOTATION_ECOLI_AMRFINDER="${ANNOTATION_ECOLI}/amrfinder/WS2762512A${sample_id}.amrfinder"
	EGGNOG_ECOLI="${ANNOTATION_ECOLI}/eggnog/WS2762512A${sample_id}.eggnog"

	echo -e "\e[31m =================== \e[0m"
	echo -e "\e[31m BAKTA: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =================== \e[0m"

	filtered_Ecoli_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"

	#functional genomic annotation, successor to prokka
	conda run -n BPannotation bakta \
	--threads ${threads} \
	--force \
	--db ${bakta_db} \
	--genus Escherichia \
	--species coli \
	--tmp-dir ${ANNOTATION_ECOLI_BAKTA} \
	--output ${ANNOTATION_ECOLI_BAKTA} \
	--prefix "WS2762512A${sample_id}.bakta" \
	${filtered_Ecoli_assembly}

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m AMRFINDER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	bakta_faa="${ANNOTATION_ECOLI_BAKTA}/WS2762512A${sample_id}.bakta.faa"
	bakta_gff="${ANNOTATION_ECOLI_BAKTA}/WS2762512A${sample_id}.bakta.gff3"

	#detect antimicrobial resistance genes, virulent factors, and stress response genes
	conda run -n ncbi amrfinder \
	--threads ${threads} \
        --plus \
        --protein "${bakta_faa}" \
        --gff "${bakta_gff}" \
        --annotation_format bakta \
        --organism Escherichia \
        --database "${amrfinder_db}" \
	--mutation_all "${ANNOTATION_ECOLI_AMRFINDER}/WS2762512A${sample_id}.mutation_all.amrfinder.tsv" \
        --output "${ANNOTATION_ECOLI_AMRFINDER}/WS2762512A${sample_id}.amrfinder.tsv"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EGGNOG-MAPPER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =========================== \e[0m"
	
	#annotate GO terms
	# Step 1: DIAMOND search only
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--dbmem \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_ECOLI}" \
	--output_dir "${EGGNOG_ECOLI}" \
	--no_annot \
	--override \
	-i "${bakta_faa}" \
	-o "WS2762512A${sample_id}.eggnog"

	# Step 2: Annotation only using hits from step 1
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_ECOLI}" \
	--output_dir "${EGGNOG_ECOLI}" \
	--annotate_hits_table "${EGGNOG_ECOLI}/WS2762512A${sample_id}.eggnog.emapper.seed_orthologs" \
	-m no_search \
	--excel \
	--tax_scope 1236 \
	--override \
	-o "WS2762512A${sample_id}.eggnog"
	
	echo -e "\e[32m Eggnog-MAPPER completed for sample WS2762512A${sample_id} \e[0m"
done

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - $(date)"
echo "========================================================"
echo "  Preprocessing       : ${PREPROCESSING_ECOLI}"
echo "  Phage assemblies    : ${ASSEMBLY_ECOLI}"
echo "  Annotation          : ${ANNOTATION_ECOLI}"
echo "  Genome completeness : ${ANNOTATION_ECOLI}/checkm"
echo "  AMR finder          : ${ANNOTATION_ECOLI}/amrfinder"
echo "  Sequence typing     : ${ANNOTATION_ECOLI}/mlst"
echo "  Full log            : ${LOG}"
echo "========================================================"

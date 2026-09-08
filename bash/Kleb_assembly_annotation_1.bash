#!/bin/bash
#scripts for de novo assembly and annotation of K. pneumoniae ILLUMINA samples
#Building No.1
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="/storage/student9/projects/bacteria_phages/Klebsiella_pneumoniae"

threads=16

#path to references and databases
truseq2="${REF_PATH}/TruSeq2-PE.fa"
truseq3="${REF_PATH}/TruSeq3-PE.fa"

bakta_db="${REF_PATH}/bakta_db/db"
amrfinder_db="${REF_PATH}/bakta_db/db/amrfinderplus-db/latest"
eggnog_db="${REF_PATH}/eggnog_db"

#create working directories

PREPROCESSING_KLEB="${WORK_PATH}/preprocessing_Kleb"
ASSEMBLY_KLEB="${WORK_PATH}/assembly_Kleb"
ANNOTATION_KLEB="${WORK_PATH}/annotation_Kleb"
ANNOTATION_KLEB_CHECKM="${ANNOTATION_KLEB}/checkm"
ANNOTATION_KLEB_MLST="${ANNOTATION_KLEB}/mlst"
QUAST_RAW_KLEB_INPUTS="${ASSEMBLY_KLEB}/quast_raw_Kleb_inputs"
CHECKM_INPUTS="${ANNOTATION_KLEB}/checkm/inputs_Kleb"

mkdir -p "${WORK_PATH}"
mkdir -p "${WORK_PATH}/preprocessing_Kleb"
mkdir -p "${WORK_PATH}/assembly_Kleb"
mkdir -p "${WORK_PATH}/annotation_Kleb"
mkdir -p "${ASSEMBLY_KLEB}/quast_raw_Kleb_inputs"
mkdir -p "${ANNOTATION_KLEB}/checkm"
mkdir -p "${ANNOTATION_KLEB}/checkm/inputs_Kleb"
mkdir -p "${ANNOTATION_KLEB}/mlst"
mkdir -p "${ANNOTATION_KLEB}/bakta"
mkdir -p "${ANNOTATION_KLEB}/amrfinder"
mkdir -p "${ANNOTATION_KLEB}/eggnog"

#create global log
LOG="${WORK_PATH}/Kleb_assembly_annotation.log"
exec > >(tee -a "${LOG}") 2>&1

echo "============================================================"
echo " K. pneumoniae assemblies and annotations - Started: $(date)"
echo "============================================================"

	#########################
	#SAMPLE DATA PREPARATION#
	#########################

	mkdir -p "${PREPROCESSING_KLEB}/fastqc_Kleb_raw"
	mkdir -p "${PREPROCESSING_KLEB}/fastqc_Kleb_pp"
		
for sample_id in "01" "02" "03" "04" "06" "07" "10" "14" "15"
do
	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: WS2762512A${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"
	
	read1="${SAMPLE_PATH}/WS2762512A${sample_id}_R1.fastq.gz"
	read2="${SAMPLE_PATH}/WS2762512A${sample_id}_R2.fastq.gz"

	#review raw fastq
	conda run -n preprocessing fastqc \
   	--threads ${threads} \
    	--outdir "${PREPROCESSING_KLEB}/fastqc_Kleb_raw" \
    	"${read1}" "${read2}"	

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m TRIMMOMATIC: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ========================= \e[0m"
	
	trim_output_1="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R1"
	trim_output_2="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R2"

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
    	--outdir "${PREPROCESSING_KLEB}/fastqc_Kleb_pp" \
	"${read1t}" "${read2t}"
done

	echo -e "\e[31m ======================================== \e[0m"
	echo -e "\e[31m MULTIQC: ALL K. PNEUMONIAE SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ======================================== \e[0m"

	#review all trimmed fastq
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_KLEB}/fastqc_Kleb_raw" \
    	--filename "${PREPROCESSING_KLEB}/multiqc_Kleb_raw"

	echo -e "\e[31m ======================================= \e[0m"
	echo -e "\e[31m MULTIQC: ALL K. PNEUMONIAE SAMPLES (PP) \e[0m"
	echo -e "\e[31m ======================================= \e[0m"

	#review all trimmed fastq
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_KLEB}/fastqc_Kleb_pp" \
    	--filename "${PREPROCESSING_KLEB}/multiqc_Kleb_pp"

	##################
	#DE NOVO ASSEMBLY#
	##################

for sample_id in "01" "02" "03" "04" "06" "07" "10" "14" "15"
do
	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SPADES: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"
	
	read1t="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R2.paired.fastq.gz"

	#de novo assembly
	conda run -n assembly spades.py \
    	-1 "${read1t}" \
    	-2 "${read2t}" \
    	--careful \
    	-t ${threads} \
    	-o "${ASSEMBLY_KLEB}/WS2762512A${sample_id}.contig.fasta"

	Kleb_assembly="${ASSEMBLY_KLEB}/WS2762512A${sample_id}.contig.fasta/contigs.fasta"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================= \e[0m"
	
	#review raw assembly statistics
	conda run -n assembly quast.py \
    	"${Kleb_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_KLEB}/WS2762512A${sample_id}.raw.quast"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SEQKIT: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"
	
	#remove contigs with less than 500bp
	conda run -n assembly seqkit seq \
    	--min-len 500 \
    	"${Kleb_assembly}" \
    	> "${ASSEMBLY_KLEB}/WS2762512A${sample_id}.contigs.filtered.fasta"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} (PP) \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	filtered_Kleb_assembly="${ASSEMBLY_KLEB}/WS2762512A${sample_id}.contigs.filtered.fasta"

	#review filtered assembly statistics
	conda run -n assembly quast.py \
    	"${filtered_Kleb_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_KLEB}/WS2762512A${sample_id}.filtered.quast"
done

	echo -e "\e[31m ====================================== \e[0m"
	echo -e "\e[31m QUAST: ALL K. PNEUMONIAE SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ====================================== \e[0m"

for sample_id in "01" "02" "03" "04" "06" "07" "10" "14" "15"
do
	cp -f "${ASSEMBLY_KLEB}/WS2762512A${sample_id}.contig.fasta/contigs.fasta" "${QUAST_RAW_KLEB_INPUTS}/WS2762512A${sample_id}.contigs.raw.fasta"
done
	#review all raw assemblies statistics
	conda run -n assembly quast.py \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_KLEB}/all_Kleb.raw.quast" \
	"${QUAST_RAW_KLEB_INPUTS}/"*.contigs.raw.fasta

	echo -e "\e[31m ===================================== \e[0m"
	echo -e "\e[31m QUAST: ALL K. PNEUMONIAE SAMPLES (PP) \e[0m"
	echo -e "\e[31m ===================================== \e[0m"

	#review all filtered assemblies statistics
	conda run -n assembly quast.py \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_KLEB}/all_Kleb.filtered.quast" \
	"${ASSEMBLY_KLEB}/"*contigs.filtered.fasta

	###########################
	#CHECK GENOME COMPLETENESS#
	###########################
	
	cp -f "${ASSEMBLY_KLEB}/"*contigs.filtered.fasta "${CHECKM_INPUTS}"

	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m CHECKM: ALL K.PNEUMONIA SAMPLES (PP) \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	#check assembly genomes completeness
	conda run -n BPstructure checkm lineage_wf \
	--threads ${threads} \
	--extension fasta \
	--tab_table \
	--reduced_tree \
	${CHECKM_INPUTS} \
	${ANNOTATION_KLEB_CHECKM}

	# CheckM quality assessment and summary
	if [ -f "${ANNOTATION_KLEB_CHECKM}/lineage.ms" ]; then
		conda run -n BPstructure checkm qa \
		"${ANNOTATION_KLEB}/checkm/lineage.ms" \
		"${ANNOTATION_KLEB_CHECKM}" \
		-o 2 \
		--tab_table \
		-f "${ANNOTATION_KLEB_CHECKM}/all_Kleb.quality.checkm.tsv"
		echo -e "\e[32m CheckM complete: ${ANNOTATION_KLEB_CHECKM}/PA445.quality.checkm.tsv \e[0m"
	else
		echo -e "\e[32m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
		echo -e "\e[32m Genome quality assessed by assembly statistics instead \e[0m"
	fi

	#############################
	#MULTI-LOCUS SEQUENCE TYPING#
	#############################

	echo -e "\e[31m =============================== \e[0m"
	echo -e "\e[31m MLST: ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m =============================== \e[0m"
	
	#serotyping strains	
	conda run -n BPtyping mlst \
	--threads ${threads} \
	--scheme klebsiella \
	--csv \
	"${ASSEMBLY_KLEB}/"*contigs.filtered.fasta \
	> "${ANNOTATION_KLEB_MLST}/WS2762512A${sample_id}.mlst.csv"

	#############
	#ANNOTATIONS#
	#############

for sample_id in "01" "02" "03" "04" "06" "07" "10" "14" "15"
do
	mkdir -p "${ANNOTATION_KLEB}/bakta/WS2762512A${sample_id}.bakta"
	mkdir -p "${ANNOTATION_KLEB}/amrfinder/WS2762512A${sample_id}.amrfinder"
	mkdir -p "${ANNOTATION_KLEB}/eggnog/WS2762512A${sample_id}.eggnog"

	ANNOTATION_KLEB_BAKTA="${ANNOTATION_KLEB}/bakta/WS2762512A${sample_id}.bakta"
	ANNOTATION_KLEB_AMRFINDER="${ANNOTATION_KLEB}/amrfinder/WS2762512A${sample_id}.amrfinder"
	EGGNOG_KLEB="${ANNOTATION_KLEB}/eggnog/WS2762512A${sample_id}.eggnog"

	echo -e "\e[31m =================== \e[0m"
	echo -e "\e[31m BAKTA: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =================== \e[0m"

	filtered_Kleb_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"

	#functional genomic annotation, successor to prokka
	conda run -n BPannotation bakta \
	--threads ${threads} \
	--force \
	--db ${bakta_db} \
	--genus Klebsiella \
	--species pneumoniae \
	--tmp-dir ${ANNOTATION_KLEB_BAKTA} \
	--output ${ANNOTATION_KLEB_BAKTA} \
	--prefix "WS2762512A${sample_id}.bakta" \
	${filtered_Kleb_assembly}

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m AMRFINDER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	bakta_faa="${ANNOTATION_KLEB_BAKTA}/WS2762512A${sample_id}.bakta.faa"
	bakta_gff="${ANNOTATION_KLEB_BAKTA}/WS2762512A${sample_id}.bakta.gff3"

	#detect antimicrobial resistance genes, virulent factors, and stress response genes
	conda run -n ncbi amrfinder \
	--threads ${threads} \
        --plus \
        --protein "${bakta_faa}" \
        --gff "${bakta_gff}" \
        --annotation_format bakta \
        --organism Klebsiella_pneumoniae \
        --database "${amrfinder_db}" \
	--mutation_all "${ANNOTATION_KLEB_AMRFINDER}/WS2762512A${sample_id}.mutation_all.amrfinder.tsv" \
        --output "${ANNOTATION_KLEB_AMRFINDER}/WS2762512A${sample_id}.amrfinder.tsv"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EGGNOG-MAPPER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =========================== \e[0m"
	
	#annotate GO terms
	# Step 1: DIAMOND search only
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--dbmem \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_KLEB}" \
	--output_dir "${EGGNOG_KLEB}" \
	--no_annot \
	--override \
	-i "${bakta_faa}" \
	-o "WS2762512A${sample_id}.eggnog"

	# Step 2: Annotation only using hits from step 1
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_KLEB}" \
	--output_dir "${EGGNOG_KLEB}" \
	--annotate_hits_table "${EGGNOG_KLEB}/WS2762512A${sample_id}.eggnog.emapper.seed_orthologs" \
	-m no_search \
	--excel \
	--tax_scope 1236 \
	--override \
	-o "WS2762512A${sample_id}.eggnog"
	
	echo -e "\e[32m Eggnog-MAPPER completed for sample WS2762512A${sample_id} \e[0m"
done

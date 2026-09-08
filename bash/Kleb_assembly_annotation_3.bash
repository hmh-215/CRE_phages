#!/bin/bash
#scripts for de novo assembly, annotation, typing, AMR and virulence profiling of K. pneumoniae ILLUMINA samples
#Building No.1 (Project 2)
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages_3"
WORK_PATH="${SAMPLE_PATH}/Klebsiella_pneumoniae"

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
AMR_KLEB="${WORK_PATH}/amr_Kleb"
ABRICATE_KLEB="${AMR_KLEB}/abricate"
PAI_KLEB="${WORK_PATH}/pai_Kleb"
KLEBORATE_KLEB="${PAI_KLEB}/kleborate"

mkdir -p "${WORK_PATH}"
mkdir -p "${PREPROCESSING_KLEB}"
mkdir -p "${ASSEMBLY_KLEB}"
mkdir -p "${ANNOTATION_KLEB}"
mkdir -p "${QUAST_RAW_KLEB_INPUTS}"
mkdir -p "${ANNOTATION_KLEB_CHECKM}"
mkdir -p "${CHECKM_INPUTS}"
mkdir -p "${ANNOTATION_KLEB_MLST}"
mkdir -p "${ANNOTATION_KLEB}/bakta"
mkdir -p "${ANNOTATION_KLEB}/amrfinder"
mkdir -p "${ANNOTATION_KLEB}/eggnog"
mkdir -p "${ABRICATE_KLEB}"
mkdir -p "${KLEBORATE_KLEB}"

sample_ids=(WS2762607A56 WS2762607A54)

#create global log
LOG="${WORK_PATH}/Kleb_assembly_annotation_3.log"
exec > >(tee -a "${LOG}") 2>&1

echo "============================================================"
echo " K. pneumoniae assemblies and annotations (Project 3) - Started: $(date)"
echo "============================================================"

	#########################
	#SAMPLE DATA PREPARATION#
	#########################

	mkdir -p "${PREPROCESSING_KLEB}/fastqc_Kleb_raw"
	mkdir -p "${PREPROCESSING_KLEB}/fastqc_Kleb_pp"

for sample_id in "${sample_ids[@]}"
do
	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: ${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	read1="${SAMPLE_PATH}/${sample_id}_R1.fastq.gz"
	read2="${SAMPLE_PATH}/${sample_id}_R2.fastq.gz"

	#review raw fastq
	conda run -n preprocessing fastqc \
   	--threads ${threads} \
    	--outdir "${PREPROCESSING_KLEB}/fastqc_Kleb_raw" \
    	"${read1}" "${read2}"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m TRIMMOMATIC: ${sample_id} \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	trim_output_1="${PREPROCESSING_KLEB}/${sample_id}.R1"
	trim_output_2="${PREPROCESSING_KLEB}/${sample_id}.R2"

	#fastq data manipulation
	conda run -n preprocessing trimmomatic PE \
	-threads ${threads} \
	-phred33 \
	 ${read1} ${read2} \
	${trim_output_1}.paired.fastq.gz ${trim_output_1}.unpaired.fastq.gz \
	${trim_output_2}.paired.fastq.gz ${trim_output_2}.unpaired.fastq.gz \
	ILLUMINACLIP:"${truseq2}":2:30:10:2:True \
	ILLUMINACLIP:"${truseq3}":2:30:10:2:True \
	LEADING:3 TRAILING:3 MINLEN:50

	read1t=${trim_output_1}.paired.fastq.gz
	read2t=${trim_output_2}.paired.fastq.gz

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m FASTQC: ${sample_id} (PP) \e[0m"
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

	#review all raw fastqc reports
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_KLEB}/fastqc_Kleb_raw" \
    	--filename "${PREPROCESSING_KLEB}/multiqc_Kleb_raw"

	echo -e "\e[31m ======================================= \e[0m"
	echo -e "\e[31m MULTIQC: ALL K. PNEUMONIAE SAMPLES (PP) \e[0m"
	echo -e "\e[31m ======================================= \e[0m"

	#review all trimmed fastqc reports
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_KLEB}/fastqc_Kleb_pp" \
    	--filename "${PREPROCESSING_KLEB}/multiqc_Kleb_pp"

	##################
	#DE NOVO ASSEMBLY#
	##################

for sample_id in "${sample_ids[@]}"
do
	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SPADES: ${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	read1t="${PREPROCESSING_KLEB}/${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_KLEB}/${sample_id}.R2.paired.fastq.gz"

	#de novo assembly
	conda run -n assembly spades.py \
    	-1 "${read1t}" \
    	-2 "${read2t}" \
	--only-assembler \
    	--isolate \
    	-t ${threads} \
    	-o "${ASSEMBLY_KLEB}/${sample_id}.contig.fasta"

	Kleb_assembly="${ASSEMBLY_KLEB}/${sample_id}.contig.fasta/contigs.fasta"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m QUAST: ${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#review raw assembly statistics
	conda run -n assembly quast.py \
    	"${Kleb_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_KLEB}/${sample_id}.raw.quast"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SEQKIT: ${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#remove contigs with less than 500bp
	conda run -n assembly seqkit seq \
    	--min-len 500 \
    	"${Kleb_assembly}" \
    	> "${ASSEMBLY_KLEB}/${sample_id}.contigs.filtered.fasta"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m QUAST: ${sample_id} (PP) \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	filtered_Kleb_assembly="${ASSEMBLY_KLEB}/${sample_id}.contigs.filtered.fasta"

	#review filtered assembly statistics
	conda run -n assembly quast.py \
    	"${filtered_Kleb_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_KLEB}/${sample_id}.filtered.quast"
done

	echo -e "\e[31m ====================================== \e[0m"
	echo -e "\e[31m QUAST: ALL K. PNEUMONIAE SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ====================================== \e[0m"

for sample_id in "${sample_ids[@]}"
do
	cp -f "${ASSEMBLY_KLEB}/${sample_id}.contig.fasta/contigs.fasta" "${QUAST_RAW_KLEB_INPUTS}/${sample_id}.contigs.raw.fasta"
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
		"${ANNOTATION_KLEB_CHECKM}/lineage.ms" \
		"${ANNOTATION_KLEB_CHECKM}" \
		-o 2 \
		--tab_table \
		-f "${ANNOTATION_KLEB_CHECKM}/all_Kleb.quality.checkm.tsv"
		echo -e "\e[32m CheckM complete: ${ANNOTATION_KLEB_CHECKM}/all_Kleb.quality.checkm.tsv \e[0m"
	else
		echo -e "\e[32m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
		echo -e "\e[32m Genome quality assessed by assembly statistics instead \e[0m"
	fi

	####################################
	#MULTI-LOCUS SEQUENCE TYPING (MLST)#
	####################################

	echo -e "\e[31m =============================== \e[0m"
	echo -e "\e[31m MLST: ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m =============================== \e[0m"

	#genotype/sequence typing
	conda run -n BPtyping mlst \
	--threads ${threads} \
	--scheme klebsiella \
	--csv \
	"${ASSEMBLY_KLEB}/"*contigs.filtered.fasta \
	> "${ANNOTATION_KLEB_MLST}/all_Kleb.mlst.csv"

	echo -e "\e[32m MLST complete: ${ANNOTATION_KLEB_MLST}/all_Kleb.mlst.csv \e[0m"

	##################################################
	#KLEBORATE: SEROTYPE (wzi/K-LOCUS) AND VIRULENCE #
	##################################################

	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m KLEBORATE: ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	#Kleborate: dedicated Klebsiella virulence and capsule/O-antigen (serotype) typing
	#Detects: yersiniabactin (ybt), colibactin (clb), aerobactin (iuc),
	#         salmochelin (iro), RmpADC (hypermucoidy), K/O locus serotype,
	#         MLST, and an AMR gene screen in one command
	conda run -n BPtyping kleborate \
	--threads ${threads} \
	--preset kpsc \
	--assemblies "${CHECKM_INPUTS}/"*.fasta \
	--outdir "${KLEBORATE_KLEB}"

	echo -e "\e[32m Kleborate complete: ${KLEBORATE_KLEB}/ \e[0m"

	#############
	#ANNOTATIONS#
	#############

for sample_id in "${sample_ids[@]}"
do
	mkdir -p "${ANNOTATION_KLEB}/bakta/${sample_id}.bakta"
	mkdir -p "${ANNOTATION_KLEB}/amrfinder/${sample_id}.amrfinder"
	mkdir -p "${ANNOTATION_KLEB}/eggnog/${sample_id}.eggnog"
	mkdir -p "${ABRICATE_KLEB}/${sample_id}.abricate"

	ANNOTATION_KLEB_BAKTA="${ANNOTATION_KLEB}/bakta/${sample_id}.bakta"
	ANNOTATION_KLEB_AMRFINDER="${ANNOTATION_KLEB}/amrfinder/${sample_id}.amrfinder"
	EGGNOG_KLEB="${ANNOTATION_KLEB}/eggnog/${sample_id}.eggnog"
	ABRICATE_OUT="${ABRICATE_KLEB}/${sample_id}.abricate"

	echo -e "\e[31m =================== \e[0m"
	echo -e "\e[31m BAKTA: ${sample_id} \e[0m"
	echo -e "\e[31m =================== \e[0m"

	filtered_Kleb_assembly="${CHECKM_INPUTS}/${sample_id}.contigs.filtered.fasta"

	#functional genomic annotation, successor to prokka
	conda run -n BPannotation bakta \
	--threads ${threads} \
	--force \
	--db ${bakta_db} \
	--genus Klebsiella \
	--species pneumoniae \
	--tmp-dir ${ANNOTATION_KLEB_BAKTA} \
	--output ${ANNOTATION_KLEB_BAKTA} \
	--prefix "${sample_id}.bakta" \
	${filtered_Kleb_assembly}

	bakta_faa="${ANNOTATION_KLEB_BAKTA}/${sample_id}.bakta.faa"
	bakta_gff="${ANNOTATION_KLEB_BAKTA}/${sample_id}.bakta.gff3"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m AMRFINDER: ${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#detect antimicrobial resistance genes, virulence factors, and stress response genes
	conda run -n ncbi amrfinder \
	--threads ${threads} \
        --plus \
        --protein "${bakta_faa}" \
        --gff "${bakta_gff}" \
        --annotation_format bakta \
        --organism Klebsiella_pneumoniae \
        --database "${amrfinder_db}" \
	--mutation_all "${ANNOTATION_KLEB_AMRFINDER}/${sample_id}.mutation_all.amrfinder.tsv" \
        --output "${ANNOTATION_KLEB_AMRFINDER}/${sample_id}.amrfinder.tsv"

	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m ABRICATE (VFDB): ${sample_id} \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	#dedicated virulence factor screen (VFDB) plus a broader AMR cross-check
	for db in vfdb card ncbi resfinder;
	do
		conda run -n BPannotation abricate \
		--db "${db}" \
		--threads ${threads} \
		--minid 80 \
		--mincov 80 \
		"${filtered_Kleb_assembly}" \
		> "${ABRICATE_OUT}/${sample_id}.abricate.${db}.tsv"

		echo -e "\e[32m ABRICATE ${db} done for ${sample_id} \e[0m"
	done

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EGGNOG-MAPPER: ${sample_id} \e[0m"
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
	-o "${sample_id}.eggnog"

	# Step 2: Annotation only using hits from step 1
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_KLEB}" \
	--output_dir "${EGGNOG_KLEB}" \
	--annotate_hits_table "${EGGNOG_KLEB}/${sample_id}.eggnog.emapper.seed_orthologs" \
	-m no_search \
	--excel \
	--tax_scope 1236 \
	--override \
	-o "${sample_id}.eggnog"

	echo -e "\e[32m Eggnog-MAPPER completed for sample ${sample_id} \e[0m"
done

	echo -e "\e[31m ============================================== \e[0m"
	echo -e "\e[31m ABRICATE: SUMMARIZE ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m ============================================== \e[0m"

	#per-database summary tables across all K. pneumoniae samples
	for db in vfdb card ncbi resfinder;
	do
		conda run -n BPannotation abricate --summary \
		"${ABRICATE_KLEB}/"*".abricate/"*".abricate.${db}.tsv" \
		> "${ABRICATE_KLEB}/all_Kleb.abricate.${db}.summary.tsv"

		echo -e "\e[32m ABRICATE summary written: all_Kleb.abricate.${db}.summary.tsv \e[0m"
	done

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - $(date)"
echo "========================================================"
echo "  Preprocessing       : ${PREPROCESSING_KLEB}"
echo "  Assemblies          : ${ASSEMBLY_KLEB}"
echo "  Annotation          : ${ANNOTATION_KLEB}"
echo "  Genome completeness : ${ANNOTATION_KLEB_CHECKM}"
echo "  Sequence typing     : ${ANNOTATION_KLEB_MLST}"
echo "  Serotype/Virulence  : ${KLEBORATE_KLEB}"
echo "  AMR (AMRFinderPlus) : ${ANNOTATION_KLEB}/amrfinder"
echo "  AMR/Virulence (ABRICATE) : ${ABRICATE_KLEB}"
echo "  Full log            : ${LOG}"
echo "========================================================"

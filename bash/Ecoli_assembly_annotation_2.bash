#!/bin/bash
#scripts for de novo assembly, annotation, typing, AMR and virulence profiling of E. coli ILLUMINA samples
#Building No.1 (Project 2)
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages_2"
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
AMR_ECOLI="${WORK_PATH}/amr_Ecoli"
ABRICATE_ECOLI="${AMR_ECOLI}/abricate"

mkdir -p "${WORK_PATH}"
mkdir -p "${PREPROCESSING_ECOLI}"
mkdir -p "${ASSEMBLY_ECOLI}"
mkdir -p "${ANNOTATION_ECOLI}"
mkdir -p "${QUAST_RAW_ECOLI_INPUTS}"
mkdir -p "${ANNOTATION_ECOLI_CHECKM}"
mkdir -p "${CHECKM_INPUTS}"
mkdir -p "${ANNOTATION_ECOLI_MLST}"
mkdir -p "${ANNOTATION_ECOLI}/bakta"
mkdir -p "${ANNOTATION_ECOLI}/amrfinder"
mkdir -p "${ANNOTATION_ECOLI}/eggnog"
mkdir -p "${ABRICATE_ECOLI}"

#create global log
LOG="${WORK_PATH}/Ecoli_assembly_annotation_2.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================"
echo " E. coli assemblies and annotations (Project 2) - Started: $(date)"
echo "======================================================"

	#########################
	#SAMPLE DATA PREPARATION#
	#########################

	mkdir -p "${PREPROCESSING_ECOLI}/fastqc_Ecoli_raw"
	mkdir -p "${PREPROCESSING_ECOLI}/fastqc_Ecoli_pp"

#E. coli samples from the ktest submission sheet (Sample # column)
sample_id=(132895-LJF30331 134194-LJF30332 144342-LJF30333 156589-LJF30334 194083-LJF30335 216381-LJF30336 216607-LJF30337 243686-LJF30338 252837-LJF30339 257625-LJF30340 319707-LJF30341 332413-LJF30342 346152-LJF30343 125919-LJF30344 131695-LJF30345 196951-LJF30346 222949-LJF30347 241208-LJF30348 251713-LJF30349 247233-LJF30350 256139-LJF30351 258569-LJF30352 EM22180-LJF30353 EM10288-LJF30354 EM21266-LJF30360 EM16897-LJF30361-W1 EM10599-LJF30364 EM10610-LJF30365)

for sample_id in "${sample_id[@]}"
do
	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: ${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	read1="${SAMPLE_PATH}/${sample_id}_L3_1.fq.gz"
	read2="${SAMPLE_PATH}/${sample_id}_L3_2.fq.gz"

	#review raw fastq
	conda run -n preprocessing fastqc \
   	--threads ${threads} \
    	--outdir "${PREPROCESSING_ECOLI}/fastqc_Ecoli_raw" \
    	"${read1}" "${read2}"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m TRIMMOMATIC: ${sample_id} \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	trim_output_1="${PREPROCESSING_ECOLI}/${sample_id}.R1"
	trim_output_2="${PREPROCESSING_ECOLI}/${sample_id}.R2"

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
    	--outdir "${PREPROCESSING_ECOLI}/fastqc_Ecoli_pp" \
	"${read1t}" "${read2t}"
done

	echo -e "\e[31m ================================== \e[0m"
	echo -e "\e[31m MULTIQC: ALL E. COLI SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ================================== \e[0m"

	#review all raw fastqc reports
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_ECOLI}/fastqc_Ecoli_raw" \
    	--filename "${PREPROCESSING_ECOLI}/multiqc_Ecoli_raw"

	echo -e "\e[31m ================================= \e[0m"
	echo -e "\e[31m MULTIQC: ALL E. COLI SAMPLES (PP) \e[0m"
	echo -e "\e[31m ================================= \e[0m"

	#review all trimmed fastqc reports
	conda run -n preprocessing multiqc \
	--force \
    	"${PREPROCESSING_ECOLI}/fastqc_Ecoli_pp" \
    	--filename "${PREPROCESSING_ECOLI}/multiqc_Ecoli_pp"

	##################
	#DE NOVO ASSEMBLY#
	##################

for sample_id in "${sample_id[@]}"
do
	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SPADES: ${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	read1t="${PREPROCESSING_ECOLI}/${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_ECOLI}/${sample_id}.R2.paired.fastq.gz"

	#de novo assembly
	conda run -n assembly spades.py \
    	-1 "${read1t}" \
    	-2 "${read2t}" \
    	--only-assembler \
    	--isolate \
    	-t ${threads} \
    	-o "${ASSEMBLY_ECOLI}/${sample_id}.contig.fasta"

	Ecoli_assembly="${ASSEMBLY_ECOLI}/${sample_id}.contig.fasta/contigs.fasta"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m QUAST: ${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#review raw assembly statistics
	conda run -n assembly quast.py \
    	"${Ecoli_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_ECOLI}/${sample_id}.raw.quast"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m SEQKIT: ${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#remove contigs with less than 500bp
	conda run -n assembly seqkit seq \
    	--min-len 500 \
    	"${Ecoli_assembly}" \
    	> "${ASSEMBLY_ECOLI}/${sample_id}.contigs.filtered.fasta"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m QUAST: ${sample_id} (PP) \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	filtered_Ecoli_assembly="${ASSEMBLY_ECOLI}/${sample_id}.contigs.filtered.fasta"

	#review filtered assembly statistics
	conda run -n assembly quast.py \
    	"${filtered_Ecoli_assembly}" \
    	--threads ${threads} \
    	--output-dir "${ASSEMBLY_ECOLI}/${sample_id}.filtered.quast"
done

	echo -e "\e[31m ================================ \e[0m"
	echo -e "\e[31m QUAST: ALL E. COLI SAMPLES (RAW) \e[0m"
	echo -e "\e[31m ================================ \e[0m"

for sample_id in "${sample_id[@]}"
do
	cp -f "${ASSEMBLY_ECOLI}/${sample_id}.contig.fasta/contigs.fasta" "${QUAST_RAW_ECOLI_INPUTS}/${sample_id}.contigs.raw.fasta"
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
		"${ANNOTATION_ECOLI_CHECKM}/lineage.ms" \
		"${ANNOTATION_ECOLI_CHECKM}" \
		-o 2 \
		--tab_table \
		-f "${ANNOTATION_ECOLI_CHECKM}/all_Ecoli.quality.checkm.tsv"
		echo -e "\e[32m CheckM complete: ${ANNOTATION_ECOLI_CHECKM}/all_Ecoli.quality.checkm.tsv \e[0m"
	else
		echo -e "\e[32m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
		echo -e "\e[32m Genome quality assessed by assembly statistics instead \e[0m"
	fi

	####################################
	#MULTI-LOCUS SEQUENCE TYPING (MLST)#
	####################################

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m MLST: ALL E. COLI SAMPLES \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#genotype/sequence typing
	conda run -n BPtyping mlst \
	--threads ${threads} \
	--scheme ecoli \
	--csv \
	"${ASSEMBLY_ECOLI}/"*contigs.filtered.fasta \
	> "${ANNOTATION_ECOLI_MLST}/all_Ecoli.mlst.csv"

	echo -e "\e[32m MLST complete: ${ANNOTATION_ECOLI_MLST}/all_Ecoli.mlst.csv \e[0m"

	#############
	#ANNOTATIONS#
	#############

for sample_id in "${sample_id[@]}"
do
	mkdir -p "${ANNOTATION_ECOLI}/bakta/${sample_id}.bakta"
	mkdir -p "${ANNOTATION_ECOLI}/amrfinder/${sample_id}.amrfinder"
	mkdir -p "${ANNOTATION_ECOLI}/eggnog/${sample_id}.eggnog"
	mkdir -p "${ABRICATE_ECOLI}/${sample_id}.abricate"

	ANNOTATION_ECOLI_BAKTA="${ANNOTATION_ECOLI}/bakta/${sample_id}.bakta"
	ANNOTATION_ECOLI_AMRFINDER="${ANNOTATION_ECOLI}/amrfinder/${sample_id}.amrfinder"
	EGGNOG_ECOLI="${ANNOTATION_ECOLI}/eggnog/${sample_id}.eggnog"
	ABRICATE_OUT="${ABRICATE_ECOLI}/${sample_id}.abricate"

	echo -e "\e[31m =================== \e[0m"
	echo -e "\e[31m BAKTA: ${sample_id} \e[0m"
	echo -e "\e[31m =================== \e[0m"

	filtered_Ecoli_assembly="${CHECKM_INPUTS}/${sample_id}.contigs.filtered.fasta"

	#functional genomic annotation, successor to prokka
	conda run -n BPannotation bakta \
	--threads ${threads} \
	--force \
	--db ${bakta_db} \
	--genus Escherichia \
	--species coli \
	--tmp-dir ${ANNOTATION_ECOLI_BAKTA} \
	--output ${ANNOTATION_ECOLI_BAKTA} \
	--prefix "${sample_id}.bakta" \
	${filtered_Ecoli_assembly}

	bakta_faa="${ANNOTATION_ECOLI_BAKTA}/${sample_id}.bakta.faa"
	bakta_gff="${ANNOTATION_ECOLI_BAKTA}/${sample_id}.bakta.gff3"

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
        --organism Escherichia \
        --database "${amrfinder_db}" \
	--mutation_all "${ANNOTATION_ECOLI_AMRFINDER}/${sample_id}.mutation_all.amrfinder.tsv" \
        --output "${ANNOTATION_ECOLI_AMRFINDER}/${sample_id}.amrfinder.tsv"

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
		"${filtered_Ecoli_assembly}" \
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
	--temp_dir "${EGGNOG_ECOLI}" \
	--output_dir "${EGGNOG_ECOLI}" \
	--no_annot \
	--override \
	-i "${bakta_faa}" \
	-o "${sample_id}.eggnog"

	# Step 2: Annotation only using hits from step 1
	conda run -n BPannotation emapper.py \
	--cpu ${threads} \
	--data_dir ${eggnog_db} \
	--temp_dir "${EGGNOG_ECOLI}" \
	--output_dir "${EGGNOG_ECOLI}" \
	--annotate_hits_table "${EGGNOG_ECOLI}/${sample_id}.eggnog.emapper.seed_orthologs" \
	-m no_search \
	--excel \
	--tax_scope 1236 \
	--override \
	-o "${sample_id}.eggnog"

	echo -e "\e[32m Eggnog-MAPPER completed for sample ${sample_id} \e[0m"
done

	echo -e "\e[31m ============================================ \e[0m"
	echo -e "\e[31m ABRICATE: SUMMARIZE ALL E. COLI SAMPLES \e[0m"
	echo -e "\e[31m ============================================ \e[0m"

	#per-database summary tables across all E. coli samples
	for db in vfdb card ncbi resfinder;
	do
		conda run -n BPannotation abricate --summary \
		"${ABRICATE_ECOLI}/"*".abricate/"*".abricate.${db}.tsv" \
		> "${ABRICATE_ECOLI}/all_Ecoli.abricate.${db}.summary.tsv"

		echo -e "\e[32m ABRICATE summary written: all_Ecoli.abricate.${db}.summary.tsv \e[0m"
	done

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - $(date)"
echo "========================================================"
echo "  Preprocessing       : ${PREPROCESSING_ECOLI}"
echo "  Assemblies          : ${ASSEMBLY_ECOLI}"
echo "  Annotation          : ${ANNOTATION_ECOLI}"
echo "  Genome completeness : ${ANNOTATION_ECOLI_CHECKM}"
echo "  Sequence typing     : ${ANNOTATION_ECOLI_MLST}"
echo "  AMR (AMRFinderPlus) : ${ANNOTATION_ECOLI}/amrfinder"
echo "  AMR/Virulence (ABRICATE) : ${ABRICATE_ECOLI}"
echo "  Full log            : ${LOG}"
echo "========================================================"
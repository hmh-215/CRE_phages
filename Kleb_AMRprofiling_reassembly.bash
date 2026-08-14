#!/bin/bash
#AMR profiling and reassembly (chromosome and plasmids) of K. pneumoniae samples
#Buiding No.2
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
SAMPLE_PATH="/storage/student9/projects/bacteria_phages"
WORK_PATH="${SAMPLE_PATH}/Klebsiella_pneumoniae"
PREPROCESSING_KLEB="${WORK_PATH}/preprocessing_Kleb"
CHECKM_INPUTS="${WORK_PATH}/annotation_Kleb/checkm/inputs_Kleb"
ANNOTATION_KLEB="${WORK_PATH}/annotation_Kleb"

threads=16
Kleb_ref="/storage/student9/references/reference_genomes/Klebsiella_pneumoniae/ncbi_dataset/data/GCF_000009885.1/GCF_000009885.1_ASM988v1_genomic.fna"
sample_id=(01 02 03 04 06 07 10 14 15)

#Paths to tools
gipsy2="/storage/student9/tools/gipsy/gipsy/gipsy2"

#Paths to databases
platon_db="${REF_PATH}/platon_db/db/"
plasmidfinder_db="${REF_PATH}/plasmidfinder_db"
card_db="${REF_PATH}/card_db"

#Create workflow directories
AMR_KLEB="${WORK_PATH}/amr_Kleb"
RGI_KLEB="${AMR_KLEB}/rgi"
ABRICATE_KLEB="${AMR_KLEB}/abricate"
MGE_KLEB="${WORK_PATH}/mge_Kleb"
ISESCAN_KLEB="${MGE_KLEB}/isescan"
INTEGRON_KLEB="${MGE_KLEB}/integron"
ICEFINDER_KLEB="${MGE_KLEB}/icefinder"
PROKKA_KLEB="${ANNOTATION_KLEB}/prokka"
PAI_KLEB="${WORK_PATH}/pai_Kleb"
GIPSY_KLEB="${PAI_KLEB}/gipsy"
KLEBORATE_KLEB="${PAI_KLEB}/kleborate"
PHISPY_KLEB="${PAI_KLEB}/phispy"
PLASMID_KLEB="${WORK_PATH}/plasmid_Kleb"
PLATON_KLEB="${PLASMID_KLEB}/platon"
PLASMIDFINDER_KLEB="${PLASMID_KLEB}/plasmidfinder"
MOBSUITE_KLEB="${PLASMID_KLEB}/mobsuite"
CHIMERA_KLEB="${WORK_PATH}/chimera_Kleb"
REASSEMBLY_KLEB="${WORK_PATH}/reassembly_Kleb"
REASSEMBLY_QUAST_INPUTS_KLEB="${REASSEMBLY_KLEB}/quast_inputs"

mkdir -p "${RGI_KLEB}"
mkdir -p "${ABRICATE_KLEB}"
mkdir -p "${ISESCAN_KLEB}"
mkdir -p "${INTEGRON_KLEB}"
mkdir -p "${ICEFINDER_KLEB}"
mkdir -p "${KLEBORATE_KLEB}"
mkdir -p "${GIPSY_KLEB}"
mkdir -p "${PHISPY_KLEB}"
mkdir -p "${PLATON_KLEB}"
mkdir -p "${PLASMIDFINDER_KLEB}"
mkdir -p "${MOBSUITE_KLEB}"
mkdir -p "${CHIMERA_KLEB}"
mkdir -p "${REASSEMBLY_KLEB}"
mkdir -p "${REASSEMBLY_QUAST_INPUTS_KLEB}"
mkdir -p "${REASSEMBLY_KLEB}/checkm_inputs"

#Global log
LOG="${WORK_PATH}/Kleb_AMRprofiling_reassembly.log"
exec > >(tee -a "${LOG}") 2>&1

echo "================================================================"
echo " K. pneumoniae AMRprofiling and reassembly analysis - Started: $(date)"
echo "================================================================"

	#########################
	#EXTENDED AMR PROFILING #
	#########################

	#RGI uses the CARD database and provides resistance mechanism,
	#gene family, and point mutation details not in AMRFinderPlus.
	#ABRICATE does a rapid screen across multiple databases for completeness.

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m RGI: LOAD CARD DATABASE \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#Load CARD reference data once before the per-sample loop
	conda run -n rgi_env rgi load \
	--card_json "${card_db}/card.json" \
	--local

	echo -e "\e[32m CARD database loaded \e[0m"

for sample_id in "${sample_id[@]}"
do
	filtered_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"
	bakta_faa="${ANNOTATION_KLEB}/bakta/WS2762512A${sample_id}.bakta/WS2762512A${sample_id}.bakta.faa"

	mkdir -p "${RGI_KLEB}/WS2762512A${sample_id}.rgi"
	mkdir -p "${ABRICATE_KLEB}/WS2762512A${sample_id}.abricate"

	RGI_OUT="${RGI_KLEB}/WS2762512A${sample_id}.rgi"
	ABRICATE_OUT="${ABRICATE_KLEB}/WS2762512A${sample_id}.abricate"

	echo -e "\e[31m ================= \e[0m"
	echo -e "\e[31m RGI: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ================= \e[0m"

	#Contig-based RGI: detects full resistance genes from assembly
	#--include_loose: capture partial hits (review manually)
	#--include_nudge: include genes near bitscore threshold
	conda run -n rgi_env rgi main \
	--input_sequence "${filtered_assembly}" \
	--output_file "${RGI_OUT}/WS2762512A${sample_id}.rgi" \
	--input_type contig \
	--alignment_tool BLAST \
	--include_loose \
	--include_nudge \
	--num_threads ${threads} \
	--clean \
	--local

	#Protein-based RGI: more sensitive for point mutations (e.g. gyrA, parC)
	#Uses Bakta protein predictions rather than 6-frame ORF translation
	conda run -n rgi_env rgi main \
	--input_sequence "${bakta_faa}" \
	--output_file "${RGI_OUT}/WS2762512A${sample_id}.rgi.protein" \
	--input_type protein \
	--alignment_tool BLAST \
	--include_loose \
	--num_threads ${threads} \
	--clean \
	--local

	echo -e "\e[32m RGI complete for WS2762512A${sample_id} \e[0m"

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m ABRICATE: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	#Screen against multiple curated databases
	for db in card ncbi resfinder argannot vfdb;
	do
		conda run -n BPannotation abricate \
		--db "${db}" \
		--threads ${threads} \
		--minid 80 \
		--mincov 80 \
		"${filtered_assembly}" \
		> "${ABRICATE_OUT}/WS2762512A${sample_id}.abricate.${db}.tsv"

		echo -e "\e[32m ABRICATE ${db} done for WS2762512A${sample_id} \e[0m"
	done
done

	echo -e "\e[31m ============================================ \e[0m"
	echo -e "\e[31m ABRICATE: SUMMARIZE ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m ============================================ \e[0m"

	#Produce per-database summary tables across all samples
	for db in card ncbi resfinder argannot vfdb;
	do
		conda run -n BPannotation abricate --summary \
		"${ABRICATE_KLEB}/"*".abricate/"*".abricate.${db}.tsv" \
		> "${ABRICATE_KLEB}/all_Kleb.abricate.${db}.summary.tsv"

		echo -e "\e[32m ABRICATE summary written: all_Kleb.abricate.${db}.summary.tsv \e[0m"
	done

	echo -e "\e[31m ========================================= \e[0m"
	echo -e "\e[31m RGI: HEATMAP OF ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m ========================================= \e[0m"

for sample_id in "${sample_id[@]}"
do
	mv "${RGI_KLEB}/WS2762512A${sample_id}.rgi/WS2762512A${sample_id}.rgi.json" "${RGI_KLEB}/WS2762512A${sample_id}rgi.json"
done

	#Aggregate RGI results into a heatmap (gene presence/absence across samples)
	conda run -n rgi_env rgi heatmap \
	--input "${RGI_KLEB}" \
	--output "${RGI_KLEB}/all_Kleb.rgi_heatmap"

	echo -e "\e[32m RGI heatmap written: ${RGI_KLEB}/all_Kleb.rgi_heatmap \e[0m"

	###############
	#MGE DETECTION#
	###############

	#ISEScan: detects IS elements (transposases, inverted terminal repeats)
	#IntegronFinder: detects class 1/2/3 integrons and gene cassettes
	#ICEfinder: web submission prepared here (no local tool available)

for sample_id in "${sample_id[@]}"
do
	filtered_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"

	mkdir -p "${ISESCAN_KLEB}/WS2762512A${sample_id}.isescan"
	mkdir -p "${INTEGRON_KLEB}/WS2762512A${sample_id}.integron"
	mkdir -p "${ICEFINDER_KLEB}/WS2762512A${sample_id}"

	ISESCAN_KLEB_SAMPLE="${ISESCAN_KLEB}/WS2762512A${sample_id}.isescan"
	INTEGRON_KLEB_SAMPLE="${INTEGRON_KLEB}/WS2762512A${sample_id}.integron"

	echo -e "\e[31m ===================== \e[0m"
	echo -e "\e[31m ISESCAN: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ===================== \e[0m"

	#Detect IS elements across the full assembly
	#Output includes: IS family, IS group, copy number, coordinates
	conda run -n recombination isescan.py \
	--nthread ${threads} \
	--seqfile "${filtered_assembly}" \
	--output "${ISESCAN_KLEB_SAMPLE}"

	echo -e "\e[32m ISEScan complete for WS2762512A${sample_id} \e[0m"

	echo -e "\e[31m ============================= \e[0m"
	echo -e "\e[31m INTEGRON FINDER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ============================= \e[0m"

	#Detect integrons and associated gene cassettes
	#--gbk: output GenBank file (compatible with JBrowse2 / Bakta)
	#--pdf: output visual map per integron
	#--circ: treat assembly contigs as circular (appropriate for complete/near-complete)
	conda run -n recombination integron_finder \
	--gbk --pdf \
	--circ \
	--outdir "${INTEGRON_KLEB_SAMPLE}" \
	"${filtered_assembly}"

	#Rename IntegronFinder outputs to consistent prefix
	INTEGRON_RESULTS="${INTEGRON_KLEB_SAMPLE}/Results_Integron_Finder_WS2762512A${sample_id}.contigs.filtered"

	if [ -d "${INTEGRON_RESULTS}" ]; then
		for contig_gbk in "${INTEGRON_RESULTS}"/contig_*.gbk;
		do
			[ -f "${contig_gbk}" ] || continue
			contig_name=$(basename "${contig_gbk}" .gbk)
			mv "${contig_gbk}" \
			"${INTEGRON_RESULTS}/WS2762512A${sample_id}.${contig_name}.integron.gbk"
		done

		pdf_count=0
		for pdf_file in "${INTEGRON_RESULTS}"/contig_*_*.pdf;
		do
			[ -f "${pdf_file}" ] || continue
			integron_num=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f1 | rev)
			contig_id=$(basename "${pdf_file}" .pdf | rev | cut -d'_' -f2- | rev | sed 's/^contig_//')
			new_name="WS2762512A${sample_id}.contig${contig_id}.integron_${integron_num}.pdf"
			mv "${pdf_file}" "${INTEGRON_RESULTS}/${new_name}"
			pdf_count=$(( pdf_count + 1 ))
		done
		echo -e "\e[32m IntegronFinder PDFs renamed: ${pdf_count} \e[0m"
	else
		echo -e "\e[32m WARNING: IntegronFinder results dir not found for WS2762512A${sample_id} \e[0m"
		echo -e "\e[32m Expected: ${INTEGRON_RESULTS} \e[0m"
	fi

	echo -e "\e[32m IntegronFinder complete for WS2762512A${sample_id} \e[0m"

	#Prepare ICEfinder submission file
	#ICEfinder has no local tool - submit via web:
	#https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html
	cp "${filtered_assembly}" \
	"${ICEFINDER_KLEB}/WS2762512A${sample_id}/WS2762512A${sample_id}_for_ICEfinder.fasta"

	echo -e "\e[32m ICEfinder submission file ready: \e[0m"
	echo -e "\e[32m   ${ICEFINDER_KLEB}/WS2762512A${sample_id}/WS2762512A${sample_id}_for_ICEfinder.fasta \e[0m"
	echo -e "\e[32m   Upload at: https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html \e[0m"
done

	###############
	#RE_ANNOTATION#
	###############

for sample_id in "${sample_id[@]}"
do
	filtered_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"
	renamed_assembly="${PROKKA_KLEB}/prokka_inputs/WS2762512A${sample_id}.renamed.fasta"	

	mkdir -p "${PROKKA_KLEB}/WS2762512A${sample_id}.prokka"
	mkdir -p "${PROKKA_KLEB}/prokka_inputs/WS2762512A${sample_id}.prokka"

	#rename contig headers
	awk '/^>/ {
		counter++
		print ">contig_" counter
		next
	}
	{ print }' "${filtered_assembly}" > "${renamed_assembly}"

	echo -e "\e[32m Renamed headers for WS2762512A${sample_id}: \e[0m"
	grep ">" "${renamed_assembly}" | head -5

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m PROKKA: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#Annotation of the K.pneumoniae strain
	conda run -n BPannotation prokka \
	--force \
	--cpus ${threads} \
	--genus Klebsiella \
	--species pneumoniae \
	--prefix WS2762512A${sample_id}.prokka \
	--outdir "${PROKKA_KLEB}/WS2762512A${sample_id}.prokka" \
	"${renamed_assembly}"

	echo -e "\e[32m Prokka complete for WS2762512A${sample_id} \e[0m"
done

	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m PROKKA: REFERENCE NTUH-K2004 \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	#Annotation of the K.pneumoniae strain
	conda run -n BPannotation prokka \
	--force \
	--cpus ${threads} \
	--genus Klebsiella \
	--species pneumoniae \
	--prefix NTUH-K2004.prokka \
	--outdir "${PROKKA_KLEB}/NTUH-K2004.prokka" \
	"${Kleb_ref}"

	echo -e "\e[32m Prokka complete for NTUH-K2004 \e[0m"


	###################################
	#PATHOGENICITY ISLANDS / VIRULENCE#
	###################################

	#GIPSy: predictions of genomic islands
	#Kleborate: Klebsiella-specific virulence loci (yersiniabactin, colibactin,
	#           aerobactin, salmochelin, RmpADC, wzi capsule typing)
	#PhiSpy: prophage detection within annotated genomes

for sample_id in "${sample_id[@]}"
do
	filtered_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"
	prokka_gbk="${ANNOTATION_KLEB}/prokka/WS2762512A${sample_id}.prokka/WS2762512A${sample_id}.prokka.gbk"
	Kleb_ref_gbk="${ANNOTATION_KLEB}/prokka/NTUH-K2004.prokka/NTUH-K2004.prokka.gbk"

	mkdir -p "${GIPSY_KLEB}/WS2762512A${sample_id}.gipsy2"
	mkdir -p "${PHISPY_KLEB}/WS2762512A${sample_id}.phispy"

	PHISPY_OUT="${PHISPY_KLEB}/WS2762512A${sample_id}.phispy"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m GIPSY2: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#export path to lib LD_LIBRARY_PATH propagation
	export LD_LIBRARY_PATH="/storage/student9/miniconda3/envs/gipsy_env/lib:${LD_LIBRARY_PATH:-}"	

	conda run -n gipsy_env ${gipsy2} \
	-q "${prokka_gbk}" \
	-s ${Kleb_ref_gbk} \
	-o "${GIPSY_KLEB}/WS2762512A${sample_id}.gipsy2" \
	-res -vir -met \
	-k fisher \
	--force

	echo -e "\e[32m GIPSY2 complete for WS2762512A${sample_id} \e[0m"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m PHISPY: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#Detect prophage sequences integrated in the chromosome
	#--output_choice 7: all output formats (coordinates, fasta, GenBank, PHASTER-style)
	#Requires GenBank input from Bakta
	conda run -n recombination PhiSpy.py \
	"${prokka_gbk}" \
	-o "${PHISPY_KLEB}/WS2762512A${sample_id}.phispy" \
	--output_choice 7 \
	--threads ${threads}

	echo -e "\e[32m PhiSpy complete for WS2762512A${sample_id} \e[0m"
done

	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m KLEBORATE: ALL K. PNEUMONIAE SAMPLES \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	#Kleborate: dedicated Klebsiella virulence and resistance typing
	#Detects: yersiniabactin (ybt), colibactin (clb), aerobactin (iuc),
	#         salmochelin (iro), RmpADC (hypermucoidy), wzi capsule type,
	#         MLST, AMR gene screen in one command
	#Running on all assemblies at once for a single summary table

	conda run -n BPtyping kleborate \
	--threads ${threads} \
	--preset kpsc \
	--assemblies "${CHECKM_INPUTS}/"*.fasta \
	--outdir "${KLEBORATE_KLEB}" 

	echo -e "\e[32m Kleborate complete: ${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv \e[0m"

	####################
	#PLASMID ASSIGNMENT#
	####################

	#Platon: contig-level classification (chromosomal vs plasmid)
	#        based on plasmid-specific marker proteins + contig properties
	#PlasmidFinder: replicon-based plasmid typing (inc group identification)
	#MOB-suite: clusters contigs into complete plasmid reconstructions
	#           based on relaxase, mate-pair formation, replicon genes

for sample_id in "${sample_id[@]}"
do
	filtered_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"

	mkdir -p "${PLATON_KLEB}/WS2762512A${sample_id}.platon"
	mkdir -p "${PLASMIDFINDER_KLEB}/WS2762512A${sample_id}.plasmidfinder"
	mkdir -p "${MOBSUITE_KLEB}/WS2762512A${sample_id}.mobsuite"

	PLATON_OUT="${PLATON_KLEB}/WS2762512A${sample_id}.platon"
	PLASMIDFINDER_OUT="${PLASMIDFINDER_KLEB}/WS2762512A${sample_id}.plasmidfinder"
	MOBSUITE_OUT="${MOBSUITE_KLEB}/WS2762512A${sample_id}.mobsuite"

	echo -e "\e[31m ==================== \e[0m"
	echo -e "\e[31m PLATON: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ==================== \e[0m"

	#Classify every contig as chromosome or plasmid
	#--mode sensitivity: include ambiguous contigs for manual review
	#Outputs: .tsv (per-contig classification), .chromosome.fasta, .plasmid.fasta
	conda run -n plasmid platon \
	--db "${platon_db}" \
	--output "${PLATON_OUT}" \
	--prefix "WS2762512A${sample_id}" \
	--mode sensitivity \
	--threads ${threads} \
	"${filtered_assembly}"

	echo -e "\e[32m Platon complete for WS2762512A${sample_id} \e[0m"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m PLASMIDFINDER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	#Identify plasmid replicons (Inc groups) - tells you which incompatibility group each plasmid belongs to
	#-x: extend search with relaxed thresholds
	conda run -n plasmid plasmidfinder.py \
	-i "${filtered_assembly}" \
	-o "${PLASMIDFINDER_OUT}" \
	-p "${plasmidfinder_db}" \
	-l 0.60 \
	-t 0.80 \
	-x

	echo -e "\e[32m PlasmidFinder complete for WS2762512A${sample_id} \e[0m"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m MOB-SUITE: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#Reconstruct complete plasmids by clustering contigs that share relaxase, mate-pair formation (mpf), and replicon markers
	#mob_recon does both typing and contig clustering in one step
	conda run -n plasmid mob_recon \
	--infile "${filtered_assembly}" \
	--outdir "${MOBSUITE_OUT}" \
	--num_threads ${threads} \
	--force

	echo -e "\e[32m MOB-suite complete for WS2762512A${sample_id} \e[0m"
	echo -e "\e[32m  chromosome.fasta     - chromosomal contigs \e[0m"
	echo -e "\e[32m  plasmid_*.fasta      - individual reconstructed plasmids \e[0m"
	echo -e "\e[32m  mobtyper_results.txt - plasmid mobility and typing \e[0m"
done

	###################
	#CHIMERA DETECTION#
	###################

	#Map original trimmed reads back to each assembly and
	#calculate per-base coverage. Regions with:
	#  (a) abrupt coverage drops to ~0 = likely misassembly breakpoint
	#  (b) very high coverage spikes = collapsed repeat / chimeric join
	#Flagstat confirms overall alignment rate (expect >95% for good assembly)

for sample_id in "${sample_id[@]}"
do
	filtered_assembly="${CHECKM_INPUTS}/WS2762512A${sample_id}.contigs.filtered.fasta"
	read1t="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R2.paired.fastq.gz"

	mkdir -p "${CHIMERA_KLEB}/WS2762512A${sample_id}.chimera"
	CHIMERA_OUT="${CHIMERA_KLEB}/WS2762512A${sample_id}.chimera"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m BWA INDEX: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#Index the filtered assembly for read mapping
	conda run -n mapping bwa index "${filtered_assembly}"

	echo -e "\e[31m ===================== \e[0m"
	echo -e "\e[31m BWA MEM: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ===================== \e[0m"

	#Map trimmed reads back to assembly; sort and index BAM
	conda run -n mapping bwa mem \
	-t ${threads} \
	"${filtered_assembly}" \
	"${read1t}" "${read2t}" \
	-o "${CHIMERA_OUT}/WS2762512A${sample_id}.tmp.sam" \
	2>> "${CHIMERA_OUT}/WS2762512A${sample_id}.bwa.log"

	echo -e "\e[31m ====================== \e[0m"
	echo -e "\e[31m SAMTOOLS: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ====================== \e[0m"

	conda run -n mapping samtools sort \
	-@ ${threads} \
	-o "${CHIMERA_OUT}/WS2762512A${sample_id}.mapped.bam" \
	"${CHIMERA_OUT}/WS2762512A${sample_id}.tmp.sam"	

	conda run -n mapping samtools index \
	"${CHIMERA_OUT}/WS2762512A${sample_id}.mapped.bam"

	#Per-contig summary statistics: mean depth, breadth of coverage
	conda run -n mapping samtools coverage \
	"${CHIMERA_OUT}/WS2762512A${sample_id}.mapped.bam" \
	> "${CHIMERA_OUT}/WS2762512A${sample_id}.coverage_summary.tsv"

	#Per-base depth across every position - required to spot
	#abrupt coverage drops indicative of chimeric joins
	conda run -n mapping samtools depth \
	-a \
	"${CHIMERA_OUT}/WS2762512A${sample_id}.mapped.bam" \
	> "${CHIMERA_OUT}/WS2762512A${sample_id}.depth_per_base.tsv"

	#Alignment summary: overall mapping rate, paired rate, chimeric pairs
	conda run -n mapping samtools flagstat \
	"${CHIMERA_OUT}/WS2762512A${sample_id}.mapped.bam" \
	> "${CHIMERA_OUT}/WS2762512A${sample_id}.flagstat.txt"

	echo -e "\e[32m Coverage stats written for WS2762512A${sample_id} \e[0m"
	echo -e "\e[32m  coverage_summary.tsv - per-contig mean depth and breadth \e[0m"
	echo -e "\e[32m  depth_per_base.tsv   - full per-base depth (load in Python/R for plotting) \e[0m"
	echo -e "\e[32m  flagstat.txt         - overall alignment rate (expect >95%) \e[0m"
done

	#Summarise flagstat mapping rates across all samples
	echo -e "\e[32m Mapping rate summary (all samples): \e[0m"

	for sample_id in "01" "02" "03" "04" "06" "07" "10" "14" "15";
	do
		rate=$(grep "mapped (" \
		"${CHIMERA_KLEB}/WS2762512A${sample_id}.chimera/WS2762512A${sample_id}.flagstat.txt" \
		| head -1 | awk '{print $5}' | tr -d '()')
		echo "  WS2762512A${sample_id}: ${rate}"
	done

	#########################
	#RE-ASSEMBLY (UNICYCLER)#
	#########################

	#Unicycler resolves circular replicons from the SPAdes assembly graph.
	#Unlike raw SPAdes, it:
	#  - Identifies and closes circular chromosomes and plasmids
	#  - Removes spurious low-depth contigs automatically
	#  - Labels output contigs as circular or linear, and with depth
	#assembly.gfa can be opened in Bandage to visualise the graph and
	#manually confirm circular topology of plasmids.

for sample_id in "${sample_id[@]}"
do
	read1t="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R1.paired.fastq.gz"
	read2t="${PREPROCESSING_KLEB}/WS2762512A${sample_id}.R2.paired.fastq.gz"

	mkdir -p "${REASSEMBLY_KLEB}/WS2762512A${sample_id}.unicycler"
	UNICYCLER_OUT="${REASSEMBLY_KLEB}/WS2762512A${sample_id}.unicycler"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m UNICYCLER: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#--mode conservative: reduces misassemblies at the cost of a few
	#more unresolved contigs -  recommended for clinical isolates
	#Output: assembly.fasta, assembly.gfa (graph), unicycler.log
	conda run -n assembly unicycler \
	-1 "${read1t}" \
	-2 "${read2t}" \
	--mode conservative \
	--threads ${threads} \
	--out "${UNICYCLER_OUT}"

	echo -e "\e[32m Unicycler complete for WS2762512A${sample_id} \e[0m"
	echo -e "\e[32m  assembly.fasta -  final contigs (header indicates circular/linear and depth) \e[0m"
	echo -e "\e[32m  assembly.gfa   -  graph for Bandage visualisation \e[0m"

	echo -e "\e[31m =================== \e[0m"
	echo -e "\e[31m QUAST: WS2762512A${sample_id} \e[0m"
	echo -e "\e[31m =================== \e[0m"

	#Assembly statistics for Unicycler output -  compare N50 and contig
	#count against the original SPAdes assembly (Building No.1)
	conda run -n assembly quast.py \
	"${UNICYCLER_OUT}/assembly.fasta" \
	--threads ${threads} \
	--output-dir "${UNICYCLER_OUT}/quast"

	echo -e "\e[32m QUAST complete for WS2762512A${sample_id} \e[0m"
done

for sample_id in "${sample_id[@]}"
do
	cp -f "${REASSEMBLY_KLEB}/WS2762512A${sample_id}.unicycler/assembly.fasta" \
	"${REASSEMBLY_QUAST_INPUTS_KLEB}/WS2762512A${sample_id}.unicycler.fasta"
done

	echo -e "\e[31m ========================================= \e[0m"
	echo -e "\e[31m QUAST: ALL UNICYCLER ASSEMBLIES (SUMMARY) \e[0m"
	echo -e "\e[31m ========================================= \e[0m"

	conda run -n assembly quast.py \
	--threads ${threads} \
	--output-dir "${REASSEMBLY_KLEB}/all_Kleb.unicycler.quast" \
	"${REASSEMBLY_QUAST_INPUTS_KLEB}/"*.fasta

	echo -e "\e[32m Multi-sample QUAST summary: ${REASSEMBLY_KLEB}/all_Kleb.unicycler.quast \e[0m"

	echo -e "\e[31m ================================================== \e[0m"
	echo -e "\e[31m CHECKM: ALL UNICYCLER ASSEMBLIES (COMPLETENESS QC) \e[0m"
	echo -e "\e[31m ================================================== \e[0m"

	#Copy Unicycler assemblies to a CheckM input directory

	for sample_id in "01" "02" "03" "04" "06" "07" "10" "14" "15";
	do
		cp "${REASSEMBLY_KLEB}/WS2762512A${sample_id}.unicycler/assembly.fasta" \
		"${REASSEMBLY_KLEB}/checkm_inputs/WS2762512A${sample_id}.unicycler.fasta"
	done

	#Re-check genome completeness after Unicycler assembly
	#Compare output against Building No.1 CheckM results
	conda run -n BPstructure checkm lineage_wf \
	-t ${threads} \
	--reduced_tree \
	--pplacer_threads 1 \
	-x fasta \
	--tab_table \
	"${REASSEMBLY_KLEB}/checkm_inputs" \
	"${REASSEMBLY_KLEB}/checkm"

	#CheckM quality assessment and summary
	if [ -f "${REASSEMBLY_KLEB}/checkm/lineage.ms" ]; then
		conda run -n BPstructure checkm qa \
		"${REASSEMBLY_KLEB}/checkm/lineage.ms" \
		"${REASSEMBLY_KLEB}/checkm" \
		-o 2 \
		--tab_table \
		-f "${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv"
		echo -e "\e[32m CheckM complete: ${REASSEMBLY_KLEB}/checkm/Kleb.reassembly.quality.checkm.tsv \e[0m"
	else
		echo -e "\e[32m WARNING: lineage.ms not found - skipping checkm qa \e[0m"
		echo -e "\e[32m Genome quality assessed by assembly statistics instead \e[0m"
	fi

echo "================================================================"
echo " K. pneumoniae follow-up analysis - Completed: $(date)"
echo "================================================================"
echo ""
echo " Output summary:"
echo "   AMR (RGI)        : ${RGI_KLEB}/"
echo "   AMR (ABRICATE)   : ${ABRICATE_KLEB}/"
echo "   IS elements      : ${ISESCAN_KLEB}/"
echo "   Integrons        : ${INTEGRON_KLEB}/"
echo "   ICEfinder inputs : ${ICEFINDER_KLEB}/"
echo "   Genomic islands  : ${ISLANDPATH_KLEB}/"
echo "   Virulence typing : ${KLEBORATE_KLEB}/all_Kleb.kleborate.tsv"
echo "   Prophages        : ${PHISPY_KLEB}/"
echo "   Platon           : ${PLATON_KLEB}/"
echo "   PlasmidFinder    : ${PLASMIDFINDER_KLEB}/"
echo "   MOB-suite        : ${MOBSUITE_KLEB}/"
echo "   Chimera check    : ${CHIMERA_KLEB}/"
echo "   Unicycler        : ${REASSEMBLY_KLEB}/"
echo "   Full log         : ${LOG}"
echo ""
echo " Manual steps required after this script:"
echo "   1. Submit ICEfinder FASTA files at:"
echo "      https://bioinfo-mml.sjtu.edu.cn/ICEfinder/ICEfinder.html"
echo "   2. Open assembly.gfa files in Bandage to confirm circular plasmids"
echo "   3. Review depth_per_base.tsv for coverage anomalies (chimera check)"
echo "   4. Cross-reference Platon + MOB-suite + PlasmidFinder to confirm"
echo "      plasmid identity and assign resistance genes to replicons"
echo "   5. Check Kleborate output for yersiniabactin/colibactin/aerobactin"
echo "      presence - report to supervisor alongside MLST"
echo "================================================================"

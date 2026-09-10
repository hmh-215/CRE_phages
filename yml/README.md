# Conda Environments for CRE_phages

This directory contains conda environment YAML files specifying all environments and tools required to run the bacterial and phage workflows in this repository.

## Environment Inventory

| Environment File | Environment Name | Tools & Key Packages | Purpose |
|---|---|---|---|
| [`preprocessing.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/preprocessing.yml) | `preprocessing` | FastQC (0.12.1), Trimmomatic (0.38), MultiQC (1.33) | Read QC and adapter trimming |
| [`mapping.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/mapping.yml) | `mapping` | BWA-MEM (0.7.19), SAMtools (1.23.1), Picard (3.4.0) | Bacterial host-read depletion, alignment |
| [`assembly.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/assembly.yml) | `assembly` | SPAdes (4.0.0), QUAST (5.2.0), seqkit (2.13.0), Unicycler (0.5.1) | De novo bacterial & phage genome assembly |
| [`BPannotation.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/BPannotation.yml) | `BPannotation` | Bakta (1.12.0), Prokka (1.15.6), Pharokka (1.9.1), PHANOTATE (1.6.7), eggNOG-mapper (2.1.13), ABRicate (1.4.0) | Bacterial & bacteriophage functional annotation |
| [`BPstructure.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/BPstructure.yml) | `BPstructure` | CheckM (1.2.5), CheckV (1.0.3) | Bacterial and viral assembly completeness / contamination QC |
| [`BPtyping.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/BPtyping.yml) | `BPtyping` | mlst (2.23.0), Kleborate (3.2.4), ECTyper (2.0.0) | Sequence typing (MLST), Klebsiella & E. coli serotyping |
| [`recombination.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/recombination.yml) | `recombination` | ISEScan (1.7.3), IntegronFinder (2.0.6), EMBOSS (6.6.0.0), PhiSpy (5.0.10) | IS elements, integrons, repeats, prophages |
| [`ncbi.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/ncbi.yml) | `ncbi` | BLAST (2.12.0), NCBI datasets CLI (18.26.0), AMRFinderPlus (4.2.7), Python helper scripts | BLAST searching, AMR detection, NCBI tools |
| [`rgi_env.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/rgi_env.yml) | `rgi_env` | RGI (6.0.5) + CARD database | Comprehensive AMR gene screening |
| [`plasmid.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/plasmid.yml) | `plasmid` | Platon (1.7), PlasmidFinder (2.1.6), MOB-suite (3.1.9) | Plasmid identification, reconstruction, and typing |
| [`phageterm_env.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/phageterm_env.yml) | `phageterm_env` | PhageTerm (4.1) | Bacteriophage termini determination |
| [`gipsy_env.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/gipsy_env.yml) | `gipsy_env` | Gipsy2 | Pathogenicity island prediction |
| [`clermontyping.yml`](file:///E:/huong/Projects/github/CRE_phages/yml/clermontyping.yml) | `clermontyping` | Clermontyping (24.02) | E. coli phylogroup assignment |

## Setup Instructions

To create an environment:
```bash
conda env create -f yml/<environment_file>.yml
```

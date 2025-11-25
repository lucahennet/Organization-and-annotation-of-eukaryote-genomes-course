# Genome annotation and comparative genomics of *Arabidopsis thaliana* accession Elh-2

## Description of the projet

This repository contains the scripts and documentation resulting from the "Organization and annotation of eukaryote genomes" course (2025) given at the University of Fribourg, Switzerland. This project directly follows the "Genome assembly and annotation"course (2025), further information about which can be found on the following GitHub page:  
https://github.com/lucahennet/Genome-Assembly-and-Annotation-Course. 

A comprehensive pipeline was executed for the annotation and comparative genomics of the *Arabidopsis thaliana* accession **Elh-2**. The main objectives were:
- Performing a systematic annotation of Transposable Elements (TEs)
- Annotating protein-coding genes by integrating homology evidence, transcriptomic data, and *ab initio* predictions
- Conducting a pangenome analysis against several related *Arabidopsis* accessions to study gene presence/absence and syntenic relationships

## Data

The assembled genome for the Elh-2 accession was obtained from the results of the "Genome assembly and annotation" course. The Hifiasm assembly was chosen specifically to ensure comparability within the class.

The datasets used in this project come from the following studies:

Lian Q, et al. A pan-genome of 69 Arabidopsis thaliana accessions reveals a conserved genome structure throughout the global species range Nature Genetics. 2024;56:982-991. Available from: https://www.nature.com/articles/s41588-024-01715-9

Jiao WB, Schneeberger K. Chromosome-level assemblies of multiple Arabidopsis genomes reveal hotspots of rearrangements with altered evolutionary dynamics. Nature Communications. 2020;11:1–10. Available from: http://dx.doi.org/10.1038/s41467-020-14779-y

Additionally, the comprehensive reference genome of Arabidopsis thaliana TAIR10 was used.

---

## Workflow and script guide

All scripts are configured as *SLURM batch jobs* (`.sh`) are designed to run on a cluster environment, often utilising Apptainer containers for dependency management. Some R scripts can be run locally. Scripts should be run in the numerical order presented in the script folder, unless explicitly notified.

    Note that the following scripts were in general diretly provided per se by the course material. However, adustments have been done to adjust the outputs or the parameters used. 

### 1. Transposable Element (TE) annotation

This phase identifies TEs and analyses TE dynamics (age and distribution) across the genome. 

#### 1.1 TE Identification and classification
* **`01.1-run_TE_annotation.sh`** - **EDTA Pipeline**    
    Executes the full EDTA (v2.2) pipeline to identify LTR retrotransposons, TIR/DNA transposons, and Helitrons. It uses the TAIR10 CDS reference to mask coding regions, performs whole-genome TE annotation (GFF3), and generates non-redundant TE library.

#### 1.2 Clade-level analysis and refinement
* **`01.2-classify_LTR_elements.sh`** - **TEsorter on raw LTRs**  
    Runs TEsorter (v1.3.0) to classify raw, full-length LTR retrotransposons, which were identified by EDTA, into specific clades (e.g., Athila, Tork) using the REXdb plant database.

* **`02.1-plot_LTR_identity_by_clade.R`** - **LTR age dynamics**  
    Reads the EDTA GFFs and TEsorter results to plot LTR identity distributions by Copia and Gypsy clade, in order to visually distinguish recent vs. ancient insertion bursts.

* **`02.2-plot_circos_TE_density.R`** - **Genome visualisation**  
    Generates circular genome plots (Circos) showing the density of TE Superfamilies across the longest scaffolds. It is recommanded to run this script after gene annotation is completed, so that a track for gene density can be added for comparative visualisation.
    
* **`03-refine_TE_classification.sh`** - **Library refinement**      
    Extracts Copia and Gypsy sequences from the non-redundant TE library and runs TEsorter separately on these subsets for amore detailed clade resolution.

#### 1.3 TE age estimation
* **`04.1-estimate_TE_age.sh`** & **`04.2-parseRM.pl`** - **Divergence analysis**  
    Executes the custom Perl script `04.2-parseRM.pl` to calculate the corrected divergence for each TE copy. This divergence value is then used to estimate the element's insertion age.

* **`04.3-plot_div.R`** - **TE landscape plot**  
    Visualises the temporal dynamics of TE accumulation (the "TE Landscape") by plotting sequence abundance against divergence.

---

### 2. Gene annotation with MAKER

This phase predicts consensus gene models by integrating all available evidence using the MAKER pipeline.

#### 2.1 MAKER setup and execution
* **`05.1-create_fai_file.sh`** - **Indexing**  
    Uses SAMtools (v1.13) to create the FASTA index (`.fai`) required by MAKER.

* **`05.2-create_control_files.sh`** - **Configuration**  
    Generates the default MAKER (v3.01.03) control files (`maker_opts.ctl`, etc.). These templates were manually edited to correctly specify the paths to the genome, transcriptome, and protein evidence paths.

* **`05.3-run-MAKER.sh`** - **Parallel prediction**       
    Executes the MAKER pipeline using MPI for parallel processing, integrating AUGUSTUS predictions, Trinity transcripts, and protein homology (TAIR10/UniProt) to generate the initial gene models.

* **`05.4-merge_MAKER_outputs.sh`** - **Consolidation**  
    Merges per-contig output files into single genome-wide GFF3 (with and without sequences) and FASTA files (proteins and transcripts).

---

### 3. Annotation refinement and functional annotation

This phase cleans the annotation (raw MAKER output), and assigns biological functions to produce the final annotation set.

#### 3.1 ID management and functional analysis
* **`06.1-rename_MAKER_ids.sh`** - **Systematic Renaming**  
    Replaces internal MAKER IDs with clean, accession-specific identifiers (e.g., `Elh-2...`) across the GFF3, protein, and transcript FASTA files.

* **`06.2-run_InterProScan_functional_annotation.sh`** - **Domain search**  
    Runs InterProScan on the predicted protein sequence to assign functional domains (Pfam) and Gene Ontology (GO) terms.

* **`06.3-filter_refine_annotation.sh`** - **Refinement**  
    Integration of the InterProScan results into the GFF3 to calculate the Annotation Edit Distance (AED) scores. Then, filtration of the low-condifence gene models, i.e., with high AED or lakc of evidence, to produce the final gene set.

#### 3.3 Homology validation
* **`08-run_homology_validation.sh`** - **BLASTP**  
    Performs BLASTP searches against UniProt and TAIR10 databases to assign functional names based on homology

---

### 4. Quality Assessment (QC)

This phase quantifies the completeness and structural statistics of the final annotation. 

* **`07.1-evaluate_annotation_completeness.sh`** - **Isoform extraction & BUSCO**  
    Extracts the single longest protein and transcript isoform for every gene, which is then used to run BUSCO (v5.4.2) against the `brassicales_odb10` lineage dataset to assess the annotation's completeness.

* **`07.2-visualise_BUSCO_qc.sh`** - **Visualisation**      
    Generates bar plots summarising the BUSCO completeness metrics for both protein and transcript sequences

* **`07.3-generate_agat_stats.sh`** - **Structural stats**  
    Runs AGAT (v1.5.1) on the final GFF3 file to generate comprehensive structural statistics, such as gene counts, exon/intron lengths, and feature distributions.

---

### 5. Comparative genomics (Pangenome)

This phase compares the annotated genome against other accessions to study pangenome structure and synteny. 

* **`09.1-setup_genespace.sh`** - **Data prep**  
    Prepares the required input files for GENESPACE by formatring the annotations of `Elh-2`, `TAIR10`and other accessions (`Are-6`, `Ice-1`, `Taz-0`) into standardised BED and FASTA files.

* **`09.2-submit_GENESPACE_run.sh`** & **`09.3-run_GENESPACE_analysis.R`** - **Pipeline execution**  
    Launches the R script, executing the core GENESPACE pipeline. The pipeline will detect syntenic blocks (MCSanX) and construct a pangenome matrix.

* **`09.4-plot_pangenome.R`** - **Analysis & plotting**     
    Analyses the pangenome matrix to classify orthogroups and genes in core (present in all), accessory and species-specific categories. Then generates final pangenome frequency plot.

---

## Dependencies
Note that the different parameters used are in general listed and explained within the scripts.

### Core tools
| Tool | Version | 
| :--- | :--- | 
| **EDTA** | 2.2 | 
| **TEsorter** | 1.3.0 | 
| **MAKER** | 3.01.03 | 
| **BUSCO** | 5.4.2 | 
| **AGAT** | 1.5.1 | 
| **InterProScan** | 5.70-102.0 | 
| **BLAST+** | 2.15.0 | 
| **GENESPACE** | Latest | 

### R packages
* `tidyverse`, `data.table`, `cowplot`
* `circlize`, `ComplexHeatmap`
* `GENESPACE`

## Maintenance Status
This project was developed as part of an academic course and will not be under active maintenance. The code is provided as-is without guarantees of future updates or support.
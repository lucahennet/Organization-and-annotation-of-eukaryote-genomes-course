#!/bin/env bash

#SBATCH --job-name=EDTA_annotation
#SBATCH --partition=pibu_el8
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/EDTA_annotation_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/EDTA_annotation_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course
CONTAINER=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif
INPUT_FASTA=$WORKDIR/hifiasm_assembly/Elh-2.asm.bp.p_ctg.fa
OUTDIR=$WORKDIR/results/EDTA_annotation

mkdir -p $OUTDIR
cd "$OUTDIR"

# Processing of the Hifiasm-assembled genome FASTA file and use a set of known CDS (Coding Sequence) data from Arabidopsis thaliana (TAIR10) as a reference for distinguishing TEs from protein-coding genes.
apptainer run --bind /data $CONTAINER \
    EDTA.pl \
    --genome $INPUT_FASTA \
    --species others \ # general mode
    --step all \ # run all steps of the EDTA pipeline
    --sensitive 1 \ # enable sensitive mode for better detection
    --cds "/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated" \ # reference coding Sequence to filter out protein-coding genes.
    --anno 1 \ # final annotation step to generate the non-redundant TE library and annotation file
    --threads $SLURM_CPUS_PER_TASK
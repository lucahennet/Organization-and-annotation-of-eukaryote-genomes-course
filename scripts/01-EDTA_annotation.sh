#!/bin/bash

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

apptainer run --bind /data $CONTAINER \
    EDTA.pl \
    --genome $INPUT_FASTA \
    --species others \
    --step all \
    --sensitive 1 \
    --cds "/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated" \
    --anno 1 \
    --threads $SLURM_CPUS_PER_TASK
#!/usr/bin/env bash

#SBATCH --job-name=run_genespace
#SBATCH --partition=pshort_el8
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/run_genespace_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/run_genespace_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/GENESPACE
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
Rscript=/data/users/lhennet/organisation_annotation_course/scripts/09.3-run_GENESPACE_analysis.R

mkdir -p $WORKDIR
cd $WORKDIR

# Execute the GENESPACE analysis from the R script
apptainer exec \
    --bind /data \
    --bind $SCRATCH:/temp \
    $COURSEDIR/containers/genespace_latest.sif Rscript $Rscript $WORKDIR
#!/bin/env bash

#SBATCH --job-name=control_file
#SBATCH --partition=pibu_el8
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/control_file_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/control_file_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation

mkdir -p $WORKDIR
cd $WORKDIR

CONTAINER=/data/courses/assembly-annotation-course/CDS_annotation/containers/MAKER_3.01.03.sif

# Generate control files for MAKER gene annotation pipeline.
apptainer exec --bind /data $CONTAINER \
    maker -CTL
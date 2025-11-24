#!/bin/env bash

#SBATCH --job-name=InterProScan
#SBATCH --partition=pibu_el8
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/InterProScan_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/InterProScan_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/final
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
protein="assembly.all.maker.proteins.fasta.renamed.fasta"

mkdir -p $WORKDIR
cd $WORKDIR

# Run InterProScan to annotate protein sequences with functional domains.
apptainer exec \
    --bind $WORKDIR:/data \
    --bind $SCRATCH:/temp \
    --bind $COURSEDIR/data/interproscan-5.70-102.0/data:/opt/interproscan/data \
    $COURSEDIR/containers/interproscan_latest.sif \
    /opt/interproscan/interproscan.sh \
    -appl pfam --disable-precalc -f TSV \ # Use Pfam application with precalculation disabled, output in TSV format
    --goterms --iprlookup --seqtype p \ # Include GO terms and InterPro annotations for protein sequences
    -i /data/$protein \ # Input protein sequences
    -o /data/output.iprscan # Output file name
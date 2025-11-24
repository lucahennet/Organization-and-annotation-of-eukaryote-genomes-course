#!/usr/bin/env bash

#SBATCH --job-name=AGAT
#SBATCH --partition=pcoursea
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/AGAT_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/AGAT_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/final
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
OUTDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/agat_stats

mkdir -p $WORKDIR
mkdir -p $OUTDIR
cd $WORKDIR

INPUTFILE="filtered.genes.renamed.gff3" # The final, quality-filtered GFF3 file

# Generate annotation statistics using AGAT
apptainer exec --bind /data \
  $COURSEDIR/containers/agat_1.5.1--pl5321hdfd78af_0.sif \
  agat_sp_statistics.pl -i "$INPUTFILE" -o $OUTDIR/annotation.stats
  # agat_sp_statistics.pl: The AGAT script to run
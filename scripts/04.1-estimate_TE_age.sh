#!/bin/env bash

#SBATCH --job-name=TE_age_estimation
#SBATCH --partition=pibu_el8
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/TE_age_estimation_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/TE_age_estimation_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/
OUTDIR=$WORKDIR/results/TE_age_estimation
SCRIPT=/data/users/lhennet/organisation_annotation_course/scripts/04.2-parseRM.pl

mkdir -p $OUTDIR
cd $OUTDIR

module add BioPerl/1.7.8-GCCcore-10.3.0

# Estimate the age of transposable elements (TEs) in the genome assembly by parsing the RepeatMasker output file using a custom Perl script.
perl $SCRIPT -i /data/users/lhennet/organisation_annotation_course/results/EDTA_annotation/Elh-2.asm.bp.p_ctg.fa.mod.EDTA.anno/Elh-2.asm.bp.p_ctg.fa.mod.out -l 50,1 -v
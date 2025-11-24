#!/bin/env bash

#SBATCH --job-name=TEsorter
#SBATCH --partition=pibu_el8
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/TEsorter_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/TEsorter_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course
CONTAINER=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif
INPUT=$WORKDIR/results/EDTA_annotation/Elh-2.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/Elh-2.asm.bp.p_ctg.fa.mod.LTR.raw.fa
OUTDIR=$WORKDIR/results/TEsorter

mkdir -p $OUTDIR
cd "$OUTDIR"

# classify FASTA file containing raw LTR (Long Terminal Repeat) elements extracted by the previous EDTA run against the RexDB-Plant database to assign them to known TE families.
apptainer run --bind /data $CONTAINER \
    TEsorter \
    $INPUT \
    -db rexdb-plant \ 
    # classification database for plants
    -p $SLURM_CPUS_PER_TASK \
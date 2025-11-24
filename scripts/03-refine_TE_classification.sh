#!/bin/bash

#SBATCH --job-name=Refining_TE_Classification
#SBATCH --partition=pibu_el8
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/Refining_TE_Classification_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/Refining_TE_Classification_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/
OUTDIR=$WORKDIR/results/TEsorter_refined
CONTAINER=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif
genome=/data/users/lhennet/organisation_annotation_course/results/EDTA_annotation/Elh-2.asm.bp.p_ctg.fa.mod.EDTA.TElib.fa

mkdir -p $OUTDIR
cd $OUTDIR

module load SeqKit/2.6.1

# Extract Copia sequences
seqkit grep -r -p "Copia" $genome > Copia_sequences.fa
# Extract Gypsy sequences
seqkit grep -r -p "Gypsy" $genome > Gypsy_sequences.fa

# Refines the classification of Copia and Gypsy LTR retrotransposon superfamilies using the TEsorter tool.
apptainer exec --bind /data $CONTAINER \
    TEsorter \
    Copia_sequences.fa \
    -db rexdb-plant \
    -p $SLURM_CPUS_PER_TASK

apptainer exec --bind /data $CONTAINER \
    TEsorter \
    Gypsy_sequences.fa \
    -db rexdb-plant \
    -p $SLURM_CPUS_PER_TASK
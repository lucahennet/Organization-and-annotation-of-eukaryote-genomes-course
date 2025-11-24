#!/bin/env bash

#SBATCH --job-name=fai_file
#SBATCH --partition=pibu_el8
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/fai_file_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/fai_file_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation
ASSEMBLY=/data/users/lhennet/organisation_annotation_course/hifiasm_assembly/Elh-2.asm.bp.p_ctg.fa

mkdir -p $WORKDIR
cd $WORKDIR

# Create FASTA index (.fai) file for the genome assembly using SAMtools.
module load SAMtools/1.13-GCC-10.3.0
samtools faidx $ASSEMBLY # Output in the hifiasm_assembly folder
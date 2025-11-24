#!/bin/bash

#SBATCH --job-name=output_preparation
#SBATCH --partition=pibu_el8
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/output_preparation_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/output_preparation_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation

mkdir -p $WORKDIR
cd $WORKDIR

MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"

# Merge MAKER output files into consolidated GFF and FASTA files

# Merge GFF with sequences
$MAKERBIN/gff3_merge -s -d \
  Elh-2.asm.bp.p_ctg.maker.output/Elh-2.asm.bp.p_ctg_master_datastore_index.log \
  > assembly.all.maker.gff

# Merge GFF without sequences
$MAKERBIN/gff3_merge -n -s -d \
  Elh-2.asm.bp.p_ctg.maker.output/Elh-2.asm.bp.p_ctg_master_datastore_index.log \
  > assembly.all.maker.noseq.gff

# Merge fasta files
$MAKERBIN/fasta_merge -d \
  Elh-2.asm.bp.p_ctg.maker.output/Elh-2.asm.bp.p_ctg_master_datastore_index.log \
  -o assembly
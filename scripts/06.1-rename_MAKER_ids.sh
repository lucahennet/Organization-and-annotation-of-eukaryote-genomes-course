#!/bin/env bash

#SBATCH --job-name=rename_genes_transcripts
#SBATCH --partition=pibu_el8
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/rename_genes_transcripts_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/rename_genes_transcripts_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"
prefix="Elh-2"

mkdir -p $WORKDIR
cd $WORKDIR

mkdir final

# Define input and output file names
protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

cp $gff final/${gff}.renamed.gff
cp $protein final/${protein}.renamed.fasta
cp $transcript final/${transcript}.renamed.fasta

cd final

# Assign clean, consistent IDs to the gene models with MAKERʼs ID mapping tools
$MAKERBIN/maker_map_ids --prefix $prefix --justify 7 ${gff}.renamed.gff \
    > id.map

$MAKERBIN/map_gff_ids id.map ${gff}.renamed.gff
$MAKERBIN/map_fasta_ids id.map ${protein}.renamed.fasta
$MAKERBIN/map_fasta_ids id.map ${transcript}.renamed.fasta
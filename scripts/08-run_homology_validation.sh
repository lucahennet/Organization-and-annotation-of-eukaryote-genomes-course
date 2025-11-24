#!/usr/bin/env bash

#SBATCH --job-name=homology_to_functionally_validated_proteins
#SBATCH --partition=pcoursea
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/homology_to_func_validated_prots_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/homology_to_func_validated_prots_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/homology_validated_proteins
mkdir -p $WORKDIR
cd $WORKDIR

MAKERBIN="/data/courses/assembly-annotation-course/CDS_annotation/softwares/Maker_v3.01.03/src/bin/"
MAKER_PROTS="/data/users/lhennet/organisation_annotation_course/gene_annotation/final/assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
UNIPROT_DB="/data/courses/assembly-annotation-course/CDS_annotation/data/uniprot/uniprot_viridiplantae_reviewed.fa"
TAIR10_DB="/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_pep_20110103_representative_gene_model"
MAKER_GFF="/data/users/lhennet/organisation_annotation_course/gene_annotation/final/filtered.genes.renamed.gff3"

UNIPROT_OUT="blastp_uniprot.outfmt6"
UNIPROT_BEST="blastp_uniprot.besthits"

# --- SEARCH 1: UNIPROT ---

module load BLAST+/2.15.0-gompi-2021a

# Run BLASTP search against UniProt
blastp -query $MAKER_PROTS -db $UNIPROT_DB -num_threads 10 -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $UNIPROT_OUT
# -outfmt 6: Standard tabular output
# -evalue 1e-5: E-value cutoff
# -max_target_seqs 10: Keep up to 10 alignments per query

# Filter the BLAST output to keep only the single best hit per query
sort -k1,1 -k12,12g $UNIPROT_OUT | sort -u -k1,1 --merge \
    > $UNIPROT_BEST
# 1. Sort by Query ID (k1,1) then by E-value (k12,12g, ascending 'g' for scientific notation)
# 2. Use 'sort -u -k1,1 --merge' to keep only the first (best-ranked) unique entry for each Query ID

# --- INTEGRATE UNIPROT RESULTS ---

# Update protein FASTA headers with functional info from the best UniProt hits
$MAKERBIN/maker_functional_fasta $UNIPROT_DB $UNIPROT_BEST $MAKER_PROTS\
    > maker_proteins.fasta.Uniprot
$MAKERBIN/maker_functional_gff $UNIPROT_DB $UNIPROT_BEST $MAKER_GFF \
    > filtered.maker.gff3.Uniprot

# --- SEARCH 2: TAIR10 ---

TAIR10_OUT="blastp_tair10.outfmt6"
TAIR10_BEST="blastp_tair10.besthits"

# 1. Create a BLAST database from the TAIR10 FASTA file
makeblastdb -in $TAIR10_DB -dbtype prot -out tair10_db

# 2. Run BLASTP search against TAIR10
blastp -query $MAKER_PROTS -db tair10_db -num_threads 10 -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $TAIR10_OUT

# 3. Filter the TAIR10 BLAST output to keep only the single best hit per query
sort -k1,1 -k12,12g $TAIR10_OUT | sort -u -k1,1 --merge \
    > $TAIR10_BEST
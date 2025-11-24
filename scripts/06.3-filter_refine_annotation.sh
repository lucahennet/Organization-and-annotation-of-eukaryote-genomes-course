#!/bin/env bash

#SBATCH --job-name=update_filter_extract
#SBATCH --partition=pibu_el8
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/update_filter_extract_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/update_filter_extract_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/final
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
MAKERBIN="$COURSEDIR/softwares/Maker_v3.01.03/src/bin"
prefix="Elh-2"

mkdir -p $WORKDIR
cd $WORKDIR

# Define variables for input files
protein="assembly.all.maker.proteins.fasta"
transcript="assembly.all.maker.transcripts.fasta"
gff="assembly.all.maker.noseq.gff"

# 3. Update GFF with InterProScan Results
echo "Updating GFF with InterProScan results..."
# Integrates functional results (domains, GO terms) into the GFF 'attributes' column
$MAKERBIN/ipr_update_gff ${gff}.renamed.gff output.iprscan \
    > ${gff}.renamed.iprscan.gff

# 4. Calculate AED Values
echo "Calculating AED values..."
# Generates a cumulative distribution function (CDF) for Annotation Edit Distance (AED)
perl $MAKERBIN/AED_cdf_generator.pl -b 0.025 ${gff}.renamed.gff \
    > assembly.all.maker.renamed.gff.AED.txt
    # -b 0.025: Sets the bin size for the histogram/CDF calculation

# 5. Filter the GFF File for Quality
echo "Filtering GFF for quality..."
# Applies a filter to remove low-quality gene models (high AED, no evidence)
perl $MAKERBIN/quality_filter.pl -s ${gff}.renamed.iprscan.gff \
    > ${gff}_iprscan_quality_filtered.gff
    # -s prints transcripts with an AED <1 and/or Pfam domain if in gff3

# 6. Filter the GFF File for Gene Features
echo "Filtering for gene features..."
# Ensures only the core structural features are kept in the final file (removes scaffold/contig lines, etc.)
grep -P "\tgene\t|\tCDS\t|\texon\t|\tfive_prime_UTR\t|\tthree_prime_UTR\t|\tmRNA\t" ${gff}_iprscan_quality_filtered.gff \
    > filtered.genes.renamed.gff3

# Check
echo "Feature types in filtered GFF:"
cut -f3 filtered.genes.renamed.gff3 | sort | uniq

# 7. Extract mRNA Sequences and Filter FASTA Files
echo "Extracting mRNA sequences..."
module load UCSC-Utils/448-foss-2021a
module load MariaDB/10.6.4-GCC-10.3.0

# Extract the unique IDs of all accepted 'mRNA' features from the final GFF
grep -P "\tmRNA\t" filtered.genes.renamed.gff3 | awk '{print $9}' | cut -d ';' -f1 | sed 's/ID=//g' \
    > list.txt

echo "Number of mRNA features found: $(wc -l < list.txt)"

# Use faSomeRecords (UCSC utility) to extract only the sequences corresponding to the accepted IDs
echo "Filtering transcript file..."
faSomeRecords ${transcript}.renamed.fasta list.txt ${transcript}.renamed.filtered.fasta
echo "Filtering protein file..."
faSomeRecords ${protein}.renamed.fasta list.txt ${protein}.renamed.filtered.fasta

echo "All steps completed successfully!"
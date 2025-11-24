#!/usr/bin/env bash

#SBATCH --job-name=BUSCO_plot
#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/BUSCO_plot_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/BUSCO_plot_%j.err

WORKDIR="/data/users/lhennet/organisation_annotation_course/gene_annotation/BUSCO"
OUTDIR=$WORKDIR/plots

mkdir -p $WORKDIR
mkdir -p $OUTDIR
cd $WORKDIR

# Load BUSCO module
module load BUSCO/5.4.2-foss-2021a

# --- Plot 1: Proteins (longest isoforms) ---
echo "Generating plot for proteins..."
# -wd: specifies the working directory containing the BUSCO run results
generate_plot.py -wd busco_proteins_longest/

# --- Plot 2: Transcripts (longest isoforms) ---
echo "Generating plot for transcripts..."
generate_plot.py -wd busco_transcripts_longest/

# --- Combined Plot ---
echo "Generating combined plot..."
COMBINED_DIR="combined_longest_summaries"
mkdir -p $COMBINED_DIR

# Copy short summary files from the two runs into the combined directory
cp busco_proteins_longest/short_summary*.txt $COMBINED_DIR/
cp busco_transcripts_longest/short_summary*.txt $COMBINED_DIR/

# Generate the plot that combines the two summary files
generate_plot.py -wd $COMBINED_DIR/

# --- MOVE PLOTS TO OUTDIR ---
echo "Moving all generated plots to $OUTDIR"
# Plots are created in $WORKDIR. This moves the plots (*_busco_summary.*) to $OUTDIR.
mv *_busco_summary.* $OUTDIR/

echo "Plot generation complete."
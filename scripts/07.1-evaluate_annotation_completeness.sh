#!/bin/env bash

#SBATCH --job-name=BUSCO
#SBATCH --partition=pibu_el8
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/BUSCO_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/BUSCO_%j.err

WORKDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/final
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation
OUTDIR=/data/users/lhennet/organisation_annotation_course/gene_annotation/BUSCO

mkdir -p $WORKDIR
mkdir -p $OUTDIR
cd $WORKDIR

# Define input files (the final, quality-filtered, renamed sequences from the previous script)
protein="assembly.all.maker.proteins.fasta.renamed.filtered.fasta"
transcript="assembly.all.maker.transcripts.fasta.renamed.filtered.fasta"

echo "Using protein file: $protein"
echo "Using transcript file: $transcript"

# --- INPUT VALIDATION ---
echo "========================================"
echo "Validating input files..."
echo "========================================"
for file in "$protein" "$transcript"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: File $file not found!"
        exit 1
    fi
    # Count the number of sequences (lines starting with '>')
    echo "$file: $(grep -c '^>' "$file") sequences"
done

# --- FUNCTION: EXTRACT LONGEST ISOFORMS ---
# Function to identify the longest protein/transcript sequence for each gene ID.
extract_longest_isoforms() {
    local input_file="$1"
    local output_file="$2"
    
    echo "Processing $input_file -> $output_file"
    
    awk '
    BEGIN {
        delete maxlen
        delete best_header
        delete best_sequence
    }
    /^>/ {
        # Process previous sequence
        if (gene != "" && seqlen > 0) {
            if (seqlen > maxlen[gene]) {
                maxlen[gene] = seqlen
                best_header[gene] = header
                best_sequence[gene] = sequence
            }
        }
        # Start new sequence
        header = $0
        gene = $0
        sub(/^>/, "", gene)
        sub(/-R.*/, "", gene)
        sequence = ""
        seqlen = 0
        next
    }
    {
        sequence = sequence $0 "\n"
        seqlen += length($0)
    }
    END {
        # Process the last sequence
        if (gene != "" && seqlen > 0) {
            if (seqlen > maxlen[gene]) {
                maxlen[gene] = seqlen
                best_header[gene] = header
                best_sequence[gene] = sequence
            }
        }
        # Output all longest isoforms
        for (g in best_header) {
            print best_header[g]
            printf "%s", best_sequence[g]
        }
    }' "$input_file" > "$output_file"
}

# --- EXECUTE ISOFORM EXTRACTION ---
echo "========================================"
echo "Extracting longest protein isoforms..."
echo "========================================"
extract_longest_isoforms "$protein" "${protein%.fasta}.longest.fasta"

echo "========================================"
echo "Extracting longest transcript isoforms..."
echo "========================================"
extract_longest_isoforms "$transcript" "${transcript%.fasta}.longest.fasta"

echo "========================================"
echo "Verifying output files..."
echo "========================================"
for file in "${protein%.fasta}.longest.fasta" "${transcript%.fasta}.longest.fasta"; do
    if [ ! -s "$file" ]; then
        echo "ERROR: Output file $file is empty!"
        exit 1
    fi
    count=$(grep -c '^>' "$file")
    echo "$file: $count sequences extracted"
    
    # Show first few headers to verify format
    echo "First 3 headers:"
    grep '^>' "$file" | head -3
    echo "---"
done

# --- BUSCO ANALYSIS ---
module load BUSCO/5.4.2-foss-2021a

# Run BUSCO on the longest protein sequences
echo "Running BUSCO on longest proteins..."
busco -i "${protein%.fasta}.longest.fasta" \
      -l brassicales_odb10 \
      -o busco_proteins_longest \
      -m proteins \
      -c $SLURM_CPUS_PER_TASK \
      --out_path $OUTDIR

# Run BUSCO on the longest transcript sequences
echo "Running BUSCO on longest transcripts..."
busco -i "${transcript%.fasta}.longest.fasta" \
      -l brassicales_odb10 \
      -o busco_transcripts_longest \
      -m transcriptome \
      -c $SLURM_CPUS_PER_TASK \
      --out_path $OUTDIR

echo "BUSCO analysis completed successfully!"
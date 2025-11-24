#!/usr/bin/env bash

#SBATCH --job-name=setup_genespace
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --output=/data/users/lhennet/organisation_annotation_course/logs/GENESPACE_preparation_%j.out
#SBATCH --error=/data/users/lhennet/organisation_annotation_course/logs/GENESPACE_preparation_%j.err

WORK_DIR="/data/users/lhennet/organisation_annotation_course/gene_annotation"
ANNOTATION_DIR="${WORK_DIR}/final"
GENESPACE_DIR="${WORK_DIR}/GENESPACE"
COURSE_DATA="/data/courses/assembly-annotation-course/CDS_annotation/data"
LIAN_GFF="${COURSE_DATA}/Lian_et_al/gene_gff/selected"
LIAN_PROTEIN="${COURSE_DATA}/Lian_et_al/protein/selected"

MY_ACCESSION="Elh-2"
OTHER_ACCESSIONS=("Are-6" "Ice-1" "Taz-0")

mkdir -p "${GENESPACE_DIR}/peptide"
mkdir -p "${GENESPACE_DIR}/bed"

# ----------------------------------------------------
# 1. PROCESS LOCAL GENOME: Elh-2
# ----------------------------------------------------
# Find the final, filtered GFF3 and longest protein FASTA files
GFF_FILE=$(find "${ANNOTATION_DIR}" -name "*filtered*.gff3" | head -n 1)
PROTEIN_FILE=$(find "${ANNOTATION_DIR}" -name "*longest*.fasta" | head -n 1)

# Create Elh-2 .bed file (Gene coordinates)
grep -P "\tgene\t" "${GFF_FILE}" | \
awk 'BEGIN{OFS="\t"} {
    split($9, a, ";");
    split(a[1], b, "=");
    gene_id = b[2];
    gsub(/[:.-]/, "_", gene_id);
    print $1, $4-1, $5, gene_id
}' | sort -k1,1 -k2,2n > "${GENESPACE_DIR}/bed/${MY_ACCESSION}.bed"

# Create Elh-2 .fa file (Peptide sequences)
# *** FIX APPLIED HERE: sub(/\..*/, "", id) to strip all dot-separated suffixes (e.g., .t1) ***
awk '
/^>/ {
    if(seq != "") {
        print seq;
    }
    id = substr($1, 2);
    sub(/-R.*/, "", id);
    sub(/\..*/, "", id); 
    gsub(/[:.-]/, "_", id);
    print ">" id;
    seq = "";
    next;
}
{
    gsub(/[^A-Za-z*]/, "", $0);
    seq = seq $0;
}
END {
    if(seq != "") {
        print seq;
    }
}
' "${PROTEIN_FILE}" > "${GENESPACE_DIR}/peptide/${MY_ACCESSION}.fa"

# ----------------------------------------------------
# 2. PROCESS TAIR10 (Reference Genome)
# ----------------------------------------------------
# Create TAIR10 .bed file (coordinates from a pre-formatted BED file)
awk 'BEGIN{OFS="\t"} {
    gene_id = $4;
    sub(/\.[0-9]+$/, "", gene_id);
    gsub(/[:.-]/, "_", gene_id);
    print $1, $2, $3, gene_id;
}' "${COURSE_DATA}/TAIR10.bed" | sort -k1,1 -k2,2n > "${GENESPACE_DIR}/bed/TAIR10.bed"

# Create TAIR10 .fa file (Peptide sequences from a pre-formatted FASTA file)
awk '
/^>/ {
    if(seq != "") {
        print seq;
    }
    id = substr($1, 2);
    sub(/\.[0-9]+$/, "", id);
    gsub(/[:.-]/, "_", id);
    print ">" id;
    seq = "";
    next;
}
{
    gsub(/[^A-Za-z*]/, "", $0);
    seq = seq $0;
}
END {
    if(seq != "") {
        print seq;
    }
}
' "${COURSE_DATA}/TAIR10.fa" > "${GENESPACE_DIR}/peptide/TAIR10.fa"

# ----------------------------------------------------
# 3. PROCESS OTHER ACCESSIONS (Loop)
# ----------------------------------------------------
for ACC in "${OTHER_ACCESSIONS[@]}"; do
    ACC_DASH="${ACC//_/-}" # Replace underscores with dashes for file path matching
    
    # Locate the correct GFF and PROT files for the current accession
    GFF="${LIAN_GFF}/${ACC_DASH}.*.gff"
    GFF=$(ls ${GFF} 2>/dev/null | head -n 1)
    
    PROT="${LIAN_PROTEIN}/${ACC_DASH}.protein.faa"
    PROT=$(ls ${PROT} 2>/dev/null | head -n 1)
    
    # Skip if either file is missing
    if [[ -z "${GFF}" ]] || [[ -z "${PROT}" ]]; then
        continue
    fi
    
    # Create .bed file (Gene coordinates) (same logic as for Elh-2 processing)
    grep -P "\tgene\t" "${GFF}" | \
    awk 'BEGIN{OFS="\t"} {
        split($9, a, ";");
        split(a[1], b, "=");
        gene_id = b[2];
        gsub(/[:.-]/, "_", gene_id);
        print $1, $4-1, $5, gene_id
    }' | sort -k1,1 -k2,2n > "${GENESPACE_DIR}/bed/${ACC}.bed"
    
    # Create .fa file (Peptide sequences)
    # *** FIX APPLIED HERE: sub(/\..*/, "", id) to strip all dot-separated suffixes (e.g., .t1) ***
    awk '
    /^>/ {
        if(seq != "") {
            print seq;
        }
        id = substr($1, 2);
        sub(/-R.*/, "", id);
        sub(/\..*/, "", id);
        gsub(/[:.-]/, "_", id);
        print ">" id;
        seq = "";
        next;
    }
    {
        gsub(/[^A-Za-z*]/, "", $0);
        seq = seq $0;
    }
    END {
        if(seq != "") {
            print seq;
        }
    }
    ' "${PROT}" > "${GENESPACE_DIR}/peptide/${ACC}.fa"
done

echo "Done"
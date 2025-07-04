#!/bin/bash

CONTIG=$1
ONT_READS=$2
MAP_PRESET="${3:-map-ont}"
MAX_ALIGNMENTS="${4:-5}"

if [ -z "$CONTIG" ] || [ -z "$ONT_READS" ]; then
    echo "Usage: $0 <contig_fasta> <ont_reads_fastq>"
    exit 1
fi
# Check if minimap2 is installed
if ! command -v minimap2 &> /dev/null; then
    echo "minimap2 could not be found. Please install it first."
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "samtools could not be found. Please install it first."
    exit 1
fi

if ! command -v paftools.js &> /dev/null; then
    echo "paftools.js could not be found. Please install it first."
    exit 1
fi

# Check if the input files exist
if [ ! -f "$CONTIG" ]; then
    echo "Contig file $CONTIG does not exist."
    exit 1
fi
if [ ! -f "$ONT_READS" ]; then
    echo "ONT reads file $ONT_READS does not exist."
    exit 1
fi
# Run minimap2 to map ONT reads to the contig
basename_contig="${CONTIG##*/}"
basename_reads="${ONT_READS##*/}"



echo "Mapping ONT reads to contig: $basename_contig"
echo "Using ONT reads from: $ONT_READS"

minimap2 -ax "$MAP_PRESET" -N "$MAX_ALIGNMENTS" "$CONTIG" "$ONT_READS" | samtools sort -o "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.sam"

# Convert SAM to PAF

paftools.js sam2paf "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.sam" > "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.paf"

#paftools incorrectly outputs CIGAR strings so recompute them
awk 'FNR==NR && !/^@/ {a[++i]=$6; next} {$18=a[++j]}1' "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.sam" "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.paf" OFS='\t' | awk '!/^@/' | sed 's/ /\t/g' > check.paf

#rename the output PAF file
mv check.paf "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.paf"


# Convert SAM to BAM
samtools view -bS "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.sam" > "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.bam"
# Index the BAM file
samtools index "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.bam"
# Clean up the sorted SAM file
rm "${basename_reads}_mapped_to_${basename_contig%.fasta}.sorted.sam"
# Print completion message
echo "Mapping completed. Output files:"

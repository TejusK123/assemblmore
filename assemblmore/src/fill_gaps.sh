#!/bin/bash

# CONTIG=$1
# ONT_READS=$2
# MAP_PRESET="${3:-map-ont}"
# MAX_ALIGNMENTS="${4:-5}"

MAKE_BAM=true
OUTPUT_DIR="$PWD"

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            echo "Usage: $0 [OPTIONS] <contig_fasta> <ont_reads_fastq> [map_preset] [max_alignments]"
            echo "Example: $0 contig.fasta reads.fastq map-ont 5"
            exit 0
            ;;
        -o|--output|--out)
            OUTPUT_DIR="${2:-$PWD}"
            shift 2
            ;;
        -b|--no-bam|--no_bam)
            MAKE_BAM=false
            ;;

        -*|--*)
            echo "Unknown option: $1"
            echo "Usage: $0 [OPTIONS] <contig_fasta> <ont_reads_fastq> [map_preset] [max_alignments]"
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            ;;
    esac
    shift
done


CONTIG="${POSITIONAL_ARGS[0]}"
ONT_READS="${POSITIONAL_ARGS[1]}"
MAP_PRESET="${POSITIONAL_ARGS[2]:-map-ont}"
MAX_ALIGNMENTS="${POSITIONAL_ARGS[3]:-5}"

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

basename_contig="${basename_contig%.fasta}"
basename_contig="${basename_contig%.fa}"
basename_contig="${basename_contig%.fna}"

echo "Mapping ONT reads to contig: $basename_contig"
echo "Using ONT reads from: $ONT_READS"
echo "Output directory: $OUTPUT_DIR"

minimap2 -ax "$MAP_PRESET" -N "$MAX_ALIGNMENTS" "$CONTIG" "$ONT_READS" | samtools sort -o "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.sam"

# Convert SAM to PAF

paftools.js sam2paf "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.sam" > "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.paf"

#paftools incorrectly outputs CIGAR strings so recompute them
awk 'FNR==NR && !/^@/ {a[++i]=$6; next} {$18=a[++j]}1' "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.sam" "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.paf" OFS='\t' | awk '!/^@/' | sed 's/ /\t/g' > "${OUTPUT_DIR}/check.paf"

#rename the output PAF file
mv "${OUTPUT_DIR}/check.paf" "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.paf"


if [ "$MAKE_BAM" = true ]; then
    echo "Converting SAM to BAM and indexing..."
    # Convert SAM to BAM
    samtools view -bS "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.sam" > "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.bam"
    # Index the BAM file
    samtools index "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.bam"
else
    echo "Skipping BAM conversion as per user request."
fi

# Clean up the sorted SAM file
rm "${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.sam"
# Print completion message
echo "Mapping completed. Output files:"
echo "  ${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.paf"
if [ "$MAKE_BAM" = true ]; then
    echo "  ${OUTPUT_DIR}/${basename_reads}_mapped_to_${basename_contig}.sorted.bam"
fi
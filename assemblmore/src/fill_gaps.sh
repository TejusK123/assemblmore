#!/bin/bash

MAKE_BAM=true
OUTPUT_DIR="$PWD"
SCRIPT_ARGS=()

# Process arguments, only intercepting script-specific options
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output|--out)
            OUTPUT_DIR="${2:-$PWD}"
            shift 2
            ;;
        -b|--no-bam|--no_bam)
            MAKE_BAM=false
            shift
            ;;
        --fill-gaps-help)
            echo "fill_gaps.sh - Wrapper script for minimap2 with PAF output generation"
            echo ""
            echo "Usage: $0 [SCRIPT_OPTIONS] [MINIMAP2_OPTIONS] <reference> <query>"
            echo ""
            echo "Script-specific options:"
            echo "  -o, --output, --out DIR    Output directory (default: current directory)"
            echo "  -b, --no-bam, --no_bam    Skip BAM file generation"
            echo "  --fill-gaps-help           Show this help message"
            echo ""
            echo "All other arguments are passed directly to minimap2."
            echo "Note: The -a flag is automatically added if not present (required for SAM output)."
            echo "For minimap2 help, use: minimap2 --help"
            echo ""
            echo "Examples:"
            echo "  $0 -x map-ont reference.fasta reads.fastq"
            echo "  $0 -o output_dir -t 8 -N 10 ref.fa reads.fq"
            echo "  $0 --no-bam -x map-pb reference.fasta reads.fastq"
            exit 0
            ;;
        *)
            # Pass all other arguments to minimap2
            SCRIPT_ARGS+=("$1")
            shift
            ;;
    esac
done

# Check if we have at least reference and query files
if [ ${#SCRIPT_ARGS[@]} -lt 2 ]; then
    echo "Error: Insufficient arguments. Need at least reference and query files."
    echo "Usage: $0 [OPTIONS] <reference> <query>"
    echo "Use --fill-gaps-help for more information."
    exit 1
fi

# Extract reference and query from the end of arguments
# Using array length to get last two elements (compatible with older bash)
ARRAY_LENGTH=${#SCRIPT_ARGS[@]}
REFERENCE="${SCRIPT_ARGS[$((ARRAY_LENGTH-2))]}"
QUERY="${SCRIPT_ARGS[$((ARRAY_LENGTH-1))]}"

# Check if -a flag is already present in arguments
HAS_SAM_OUTPUT=false
for arg in "${SCRIPT_ARGS[@]:0:$((ARRAY_LENGTH-2))}"; do
    if [ "$arg" = "-a" ]; then
        HAS_SAM_OUTPUT=true
        break
    fi
done

# If -a is not present, we need to add it for SAM output
if [ "$HAS_SAM_OUTPUT" = false ]; then
    # Insert -a at the beginning of minimap2 arguments (before reference and query)
    MINIMAP2_ARGS=("-a" "${SCRIPT_ARGS[@]:0:$((ARRAY_LENGTH-2))}" "$REFERENCE" "$QUERY")
else
    # Use arguments as-is
    MINIMAP2_ARGS=("${SCRIPT_ARGS[@]}")
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

# Check if input files exist
if [ ! -f "$REFERENCE" ]; then
    echo "Reference file $REFERENCE does not exist."
    exit 1
fi
if [ ! -f "$QUERY" ]; then
    echo "Query file $QUERY does not exist."
    exit 1
fi

# Generate output filenames
basename_reference="${REFERENCE##*/}"
basename_query="${QUERY##*/}"

# Remove common file extensions
for ext in .fasta .fa .fna .fastq .fq; do
    basename_reference="${basename_reference%$ext}"
    basename_query="${basename_query%$ext}"
done

echo "Mapping with minimap2..."
echo "Reference: $REFERENCE"
echo "Query: $QUERY"
echo "Output directory: $OUTPUT_DIR"

# Run minimap2 with all passed arguments
echo "Running: minimap2 ${MINIMAP2_ARGS[*]} | samtools sort -o ${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.sam"
minimap2 "${MINIMAP2_ARGS[@]}" | samtools sort -o "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.sam"


# Convert SAM to PAF
paftools.js sam2paf "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.sam" > "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.paf"

# paftools incorrectly outputs CIGAR strings so recompute them
awk 'FNR==NR && !/^@/ {a[++i]=$6; next} {$18=a[++j]}1' "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.sam" "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.paf" OFS='\t' | awk '!/^@/' | sed 's/ /\t/g' > "${OUTPUT_DIR}/check.paf"

# Rename the output PAF file
mv "${OUTPUT_DIR}/check.paf" "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.paf"

if [ "$MAKE_BAM" = true ]; then
    echo "Converting SAM to BAM and indexing..."
    # Convert SAM to BAM
    samtools view -bS "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.sam" > "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.bam"
    # Index the BAM file
    samtools index "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.bam"
else
    echo "Skipping BAM conversion as per user request."
fi

# Clean up the sorted SAM file
rm "${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.sam"

# Print completion message
echo "Mapping completed. Output files:"
echo "  ${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.paf"
if [ "$MAKE_BAM" = true ]; then
    echo "  ${OUTPUT_DIR}/${basename_query}_mapped_to_${basename_reference}.sorted.bam"
fi
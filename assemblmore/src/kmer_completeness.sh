#!/bin/bash

# K-mer completeness analysis using meryl
# This script computes assembly completeness based on k-mer presence in raw reads vs assembly

# Check dependencies
if ! command -v meryl >/dev/null 2>&1; then
    echo "Error: meryl is not installed. Please install meryl to use k-mer completeness analysis."
    exit 1
fi

if ! command -v seqkit >/dev/null 2>&1; then
    echo "Error: seqkit is not installed. Please install seqkit for sequence statistics."
    exit 1
fi

if ! command -v bc >/dev/null 2>&1; then
    echo "Error: bc is not installed. Please install bc for mathematical calculations."
    exit 1
fi

if [ -z "$MERQURY" ]; then
    echo "Warning: MERQURY environment variable is not set."
    echo "Will use heuristic k-mer size calculation instead of best_k.sh"
    USE_HEURISTIC=true
else
    if [ ! -f "$MERQURY/best_k.sh" ]; then
        echo "Warning: best_k.sh not found in MERQURY directory: $MERQURY"
        echo "Will use heuristic k-mer size calculation instead"
        USE_HEURISTIC=true
    else
        USE_HEURISTIC=false
    fi
fi

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <raw_reads.fastq> <assembly.fasta>"
    echo "  raw_reads.fastq: Raw sequencing reads file"
    echo "  assembly.fasta:  Assembly file to analyze for completeness"
    exit 1
fi

RAW_READS="$1"
ASSEMBLY="$2"

# Check if the input files exist
if [ ! -f "$RAW_READS" ]; then
    echo "Error: Raw reads file '$RAW_READS' does not exist."
    exit 1
fi

if [ ! -f "$ASSEMBLY" ]; then
    echo "Error: Assembly file '$ASSEMBLY' does not exist."
    exit 1
fi

# Create a temporary directory for meryl databases
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

echo "=== K-mer Completeness Analysis ==="
echo "Raw reads: $RAW_READS"
echo "Assembly: $ASSEMBLY"
echo "Working directory: $TMP_DIR"
echo

# Step 1: Calculate optimal k-mer size using assembly length
echo "Step 1: Calculating optimal k-mer size..."
ASSEMBLY_LENGTH=$(seqkit stats -T "$ASSEMBLY" | tail -n +2 | cut -f5 | awk '{sum += $1} END {print sum}')
echo "Total assembly length: $ASSEMBLY_LENGTH bp"

# Use merqury's best_k.sh to get optimal k-mer size, or use heuristic
if [ "$USE_HEURISTIC" = "true" ]; then
    # Simple heuristic: log4(genome_size) + 10, capped between 15 and 31
    OPTIMAL_K=$(echo "scale=0; l($ASSEMBLY_LENGTH)/l(4) + 10" | bc -l | awk '{print int($1)}')
    if [ "$OPTIMAL_K" -lt 15 ]; then
        OPTIMAL_K=15
    elif [ "$OPTIMAL_K" -gt 31 ]; then
        OPTIMAL_K=31
    fi
    echo "Using heuristic k-mer size calculation"
else
    OPTIMAL_K_OUTPUT=$("$MERQURY/best_k.sh" "$ASSEMBLY_LENGTH" 2>/dev/null | tail -1)
    OPTIMAL_K=$(echo "$OPTIMAL_K_OUTPUT" | awk '{print int($1)}')  # Floor of the output
    echo "Using merqury best_k.sh calculation"
fi

if [ -z "$OPTIMAL_K" ] || [ "$OPTIMAL_K" -le 0 ]; then
    echo "Warning: Could not determine optimal k-mer size, using default k=21"
    OPTIMAL_K=21
fi

echo "Optimal k-mer size: $OPTIMAL_K"
echo

# Get base names for output files
RAW_READS_BASENAME=$(basename "$RAW_READS" | sed 's/\.[^.]*$//')
ASSEMBLY_BASENAME=$(basename "$ASSEMBLY" | sed 's/\.[^.]*$//')

# Step 2: Count k-mers in raw reads
echo "Step 2: Counting k-mers in raw reads..."
RAW_READS_MERYL="$TMP_DIR/${RAW_READS_BASENAME}.meryl"
meryl k="$OPTIMAL_K" count "$RAW_READS" output "$RAW_READS_MERYL"
if [ $? -ne 0 ]; then
    echo "Error: Failed to count k-mers in raw reads"
    exit 1
fi
echo "Raw reads k-mer database created: $RAW_READS_MERYL"
echo

# Step 3: Count k-mers in assembly
echo "Step 3: Counting k-mers in assembly..."
ASSEMBLY_MERYL="$TMP_DIR/${ASSEMBLY_BASENAME}.meryl"
meryl k="$OPTIMAL_K" count "$ASSEMBLY" output "$ASSEMBLY_MERYL"
if [ $? -ne 0 ]; then
    echo "Error: Failed to count k-mers in assembly"
    exit 1
fi
echo "Assembly k-mer database created: $ASSEMBLY_MERYL"
echo

# Step 4: Calculate difference (k-mers in reads but not in assembly)
echo "Step 4: Calculating k-mers missing from assembly..."
DIFFERENCE_MERYL="$TMP_DIR/missing_kmers.meryl"
meryl difference "$RAW_READS_MERYL" "$ASSEMBLY_MERYL" output "$DIFFERENCE_MERYL"
if [ $? -ne 0 ]; then
    echo "Error: Failed to calculate k-mer difference"
    exit 1
fi
echo "Missing k-mers database created: $DIFFERENCE_MERYL"
echo

# Step 5: Calculate intersection-sum (k-mers present in both)
echo "Step 5: Calculating k-mers present in both reads and assembly..."
INTERSECTION_MERYL="$TMP_DIR/shared_kmers.meryl"
meryl intersect-sum "$RAW_READS_MERYL" "$ASSEMBLY_MERYL" output "$INTERSECTION_MERYL"
if [ $? -ne 0 ]; then
    echo "Error: Failed to calculate k-mer intersection"
    exit 1
fi
echo "Shared k-mers database created: $INTERSECTION_MERYL"
echo

# Step 6: Calculate completeness metric
echo "Step 6: Calculating k-mer completeness..."

# Extract statistics for shared k-mers (x1)
X1=$(meryl statistics "$INTERSECTION_MERYL" | head -n4 | tail -n1 | awk '{print $2}')
if [ -z "$X1" ]; then
    echo "Error: Could not extract shared k-mer count"
    exit 1
fi

# Extract statistics for missing k-mers (x2)  
X2=$(meryl statistics "$DIFFERENCE_MERYL" | head -n4 | tail -n1 | awk '{print $2}')
if [ -z "$X2" ]; then
    echo "Error: Could not extract missing k-mer count"
    exit 1
fi

# Calculate completeness: x1 / (x1 + x2)
TOTAL_KMERS=$((X1 + X2))
if [ "$TOTAL_KMERS" -eq 0 ]; then
    echo "Error: No k-mers found in analysis"
    exit 1
fi

COMPLETENESS=$(echo "scale=6; $X1 / ($X1 + $X2)" | bc -l)
COMPLETENESS_PERCENT=$(echo "scale=2; $COMPLETENESS * 100" | bc -l)

echo "=== K-mer Completeness Results ==="
echo "K-mer size used: $OPTIMAL_K"
echo "Shared k-mers (present in both): $X1"
echo "Missing k-mers (in reads only): $X2"
echo "Total k-mers in reads: $TOTAL_KMERS"
echo "K-mer completeness: $COMPLETENESS ($COMPLETENESS_PERCENT%)"
echo

# Save results to file
RESULTS_FILE="kmer_completeness_${ASSEMBLY_BASENAME}_results.txt"
cat > "$RESULTS_FILE" << EOF
K-mer Completeness Analysis Results
===================================
Date: $(date)
Raw reads file: $RAW_READS
Assembly file: $ASSEMBLY
K-mer size: $OPTIMAL_K
Assembly length: $ASSEMBLY_LENGTH bp

Results:
--------
Shared k-mers (in both reads and assembly): $X1
Missing k-mers (in reads but not assembly): $X2
Total k-mers in reads: $TOTAL_KMERS
K-mer completeness: $COMPLETENESS ($COMPLETENESS_PERCENT%)

Interpretation:
--------------
This metric represents the fraction of k-mers from the raw reads
that are also present in the assembly. Higher values indicate
better assembly completeness.
EOF

echo "Results saved to: $RESULTS_FILE"
echo "Analysis complete!"
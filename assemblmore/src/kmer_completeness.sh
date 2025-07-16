#!/bin/bash

# K-mer completeness analysis using meryl
# This script computes assembly completeness based on k-mer presence in raw reads vs assembly

# Check dependencies

check_dependencies() {
    local dependencies=("meryl" "seqkit" "bc")
    for dep in "${dependencies[@]}"; do
        if ! command -v "$dep" >/dev/null 2>&1; then
            echo "Error: $dep is not installed. Please install $dep to use k-mer completeness analysis."
            exit 1
        fi
    done

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
}


gen_dbs() {
    local fasta_file="$1"
    local kmer_size="$2"
    local output_db="$3"

    echo "Counting k-mers in $fasta_file with k=$kmer_size..."
    
    # Check if raw database already exists
    if [ ! -d "$output_db.meryl" ]; then
        # Generate new database
        meryl k="$kmer_size" count "$fasta_file" output "${output_db}.meryl"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to count k-mers in $fasta_file"
            return 1
        fi
        echo "Raw k-mer database created: ${output_db}.meryl"
    else
        echo "Found existing k-mer database: ${output_db}.meryl"
    fi

    # Check if histogram and ploidy files exist, generate if needed
    if [ ! -f "${output_db}.hist" ]; then
        echo "Generating histogram for $output_db..."
        meryl histogram "$output_db" > "${output_db}.hist"
    fi

    if [ ! -f "${output_db}.hist.ploidy" ]; then
        echo "Generating ploidy depth analysis for $output_db..."
        java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.jar "${output_db}.hist" > "${output_db}.hist.ploidy"
    fi

    # Get filter value
    filt=`sed -n 2p "${output_db}.hist.ploidy" | awk '{print $NF}'`
    #2, default, 3, haploid, 4, diploid <- its hardcoded in merqury so whats the point of the three values????
    # Check if filtered database exists, generate if needed
    if [ ! -d "${output_db}.gt${filt}.meryl" ]; then
        echo "Generating filtered database for $output_db..."
        meryl greater-than $filt output "${output_db}.gt${filt}.meryl" "${output_db}.meryl"
        echo "Filtered k-mer database created: ${output_db}.gt${filt}.meryl"
    else
        echo "Found existing filtered database: ${output_db}.gt${filt}.meryl"
    fi
    
    return 0
}

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

check_dependencies


# Create a persistent directory for meryl databases
KMER_DB_DIR="kmer_db"
mkdir -p "$KMER_DB_DIR"

echo "=== K-mer Completeness Analysis ==="
echo "Raw reads: $RAW_READS"
echo "Assembly: $ASSEMBLY"
echo "K-mer database directory: $KMER_DB_DIR"
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
RAW_READS_MERYL="$KMER_DB_DIR/${RAW_READS_BASENAME%%.*}_k${OPTIMAL_K}"
echo "Generating k-mer database for raw reads: $RAW_READS_MERYL"
gen_dbs "$RAW_READS" "$OPTIMAL_K" "$RAW_READS_MERYL"
if [ $? -ne 0 ]; then
    echo "Error: Failed to generate k-mer database for raw reads"
    exit 1
fi

# Determine the filtered database path for reads
READS_FILT=$(sed -n 2p "${RAW_READS_MERYL%.meryl}.hist.ploidy" | awk '{print $NF}')
RAW_READS_FILTERED="$KMER_DB_DIR/${RAW_READS_BASENAME%%.*}_k${OPTIMAL_K}.gt${READS_FILT}.meryl"
echo "Using filtered reads database: $RAW_READS_FILTERED"
echo

# Step 3: Count k-mers in assembly 
echo "Step 3: Counting k-mers in assembly..."
ASSEMBLY_MERYL="$KMER_DB_DIR/${ASSEMBLY_BASENAME%%.*}_k${OPTIMAL_K}"
echo "Generating k-mer database for assembly: $ASSEMBLY_MERYL"
gen_dbs "$ASSEMBLY" "$OPTIMAL_K" "$ASSEMBLY_MERYL"
if [ $? -ne 0 ]; then
    echo "Error: Failed to generate k-mer database for assembly"
    exit 1
fi

# Determine the filtered database path for assembly
ASSEMBLY_FILT=$(sed -n 2p "${ASSEMBLY_MERYL%.meryl}.hist.ploidy" | awk '{print $NF}')
ASSEMBLY_FILTERED="$KMER_DB_DIR/${ASSEMBLY_BASENAME%%.*}_k${OPTIMAL_K}.gt${ASSEMBLY_FILT}.meryl"
echo "Using filtered assembly database: $ASSEMBLY_FILTERED"
echo

# Step 4: Calculate difference (k-mers in reads but not in assembly)
echo "Step 4: Calculating k-mers missing from assembly..."
READS_DIFF_ASSEMBLY="$KMER_DB_DIR/${RAW_READS_BASENAME}_diff_${ASSEMBLY_BASENAME}_k${OPTIMAL_K}_filtered.meryl"
if [ -d "$READS_DIFF_ASSEMBLY" ]; then
    echo "Found existing difference database: $READS_DIFF_ASSEMBLY"
    echo "Skipping difference calculation (using existing database)"
else
    meryl difference "$RAW_READS_FILTERED" "$ASSEMBLY_FILTERED" output "$READS_DIFF_ASSEMBLY"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to calculate k-mer difference"
        exit 1
    fi
    echo "Missing k-mers database created: $READS_DIFF_ASSEMBLY"
fi
echo

ASSEMBLY_DIFF_READS="$KMER_DB_DIR/${ASSEMBLY_BASENAME}_diff_${RAW_READS_BASENAME}_k${OPTIMAL_K}_filtered.meryl"
if [ -d "$ASSEMBLY_DIFF_READS" ]; then
    echo "Found existing assembly difference database: $ASSEMBLY_DIFF_READS"
    echo "Skipping assembly difference calculation (using existing database)"
else
    meryl difference "$ASSEMBLY_FILTERED" "$RAW_READS_FILTERED" output "$ASSEMBLY_DIFF_READS"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to calculate k-mers missing from reads"
        exit 1
    fi
    echo "Assembly missing k-mers database created: $ASSEMBLY_DIFF_READS"
fi
echo


# Step 5: Calculate intersection-sum (k-mers present in both)
echo "Step 5: Calculating k-mers present in both reads and assembly..."
INTERSECTION_MERYL="$KMER_DB_DIR/${RAW_READS_BASENAME}_intersect_${ASSEMBLY_BASENAME}_k${OPTIMAL_K}_filtered.meryl"
if [ -d "$INTERSECTION_MERYL" ]; then
    echo "Found existing intersection database: $INTERSECTION_MERYL"
    echo "Skipping intersection calculation (using existing database)"
else
    meryl intersect  "$ASSEMBLY_FILTERED" "$RAW_READS_FILTERED" output "$INTERSECTION_MERYL"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to calculate k-mer intersection"
        exit 1
    fi
    echo "Shared k-mers database created: $INTERSECTION_MERYL"
fi
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
X2=$(meryl statistics "$READS_DIFF_ASSEMBLY" | head -n4 | tail -n1 | awk '{print $2}')
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

X_READS=$(meryl statistics "$RAW_READS_FILTERED" | head -n4 | tail -n1 | awk '{print $2}')
if [ -z "$X_READS" ]; then
    echo "Error: Could not extract total k-mers in reads"
    exit 1
fi

COMPLETENESS=$(echo "scale=6; $X1 / ($X1 + $X2)" | bc -l)
COMPLETENESS_PERCENT=$(echo "scale=2; $COMPLETENESS * 100" | bc -l)

COMPLETENESS2=$(echo "scale=6; $X1 / $X_READS" | bc -l)
COMPLETENESS_PERCENT2=$(echo "scale=2; $COMPLETENESS2 * 100" | bc -l)

X3=$(meryl statistics "$ASSEMBLY_DIFF_READS" | head -n4 | tail -n1 | awk '{print $2}')
if [ -z "$X3" ]; then
    echo "Error: Could not extract assembly missing k-mer count"
    exit 1
fi

X4=$(meryl statistics "$ASSEMBLY_FILTERED" | head -n4 | tail -n1 | awk '{print $2}')
if [ -z "$X4" ]; then
    echo "Error: Could not extract total k-mer count in assembly"
    exit 1
fi

ERROR=$(echo "$X3 $X4" | awk -v k=$OPTIMAL_K '{print (1-(1-$1/$2)^(1/k))}')
QV=$(echo "$X3 $X4" | awk -v k=$OPTIMAL_K '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}')

echo "=== K-mer Completeness Results ==="
echo "K-mer size used: $OPTIMAL_K"
echo "Shared k-mers (present in both): $X1"
echo "Missing k-mers (in reads only): $X2"
echo "Total k-mers in reads: $TOTAL_KMERS"
echo "K-mer completeness: $COMPLETENESS ($COMPLETENESS_PERCENT%)"
echo "K-mer completeness (relative to reads): $COMPLETENESS2 ($COMPLETENESS_PERCENT2%)"
echo "Error Rate: $ERROR"
echo "Quality Value (QV): $QV"
echo

# Save results to file
RESULTS_FILE="kmer_completeness_${ASSEMBLY_BASENAME}_relative_${RAW_READS_BASENAME}.txt"
cat > "$RESULTS_FILE" << EOF
K-mer Completeness Analysis Results
===================================
Date: $(date)
Raw reads file: $RAW_READS
Assembly file: $ASSEMBLY
K-mer size: $OPTIMAL_K
Assembly length: $ASSEMBLY_LENGTH bp
K-mer database directory: $KMER_DB_DIR

Results:
--------
Shared k-mers (in both reads and assembly): $X1
Missing k-mers (in reads but not assembly): $X2
Total k-mers in reads: $TOTAL_KMERS
K-mer completeness: $COMPLETENESS ($COMPLETENESS_PERCENT%)
K-mer completeness (relative to reads): $COMPLETENESS2 ($COMPLETENESS_PERCENT2%)

Error Rate: $ERROR
Quality Value (QV): $QV

Database Files:
--------------
Raw reads k-mer database: $RAW_READS_MERYL
Raw reads filtered database: $RAW_READS_FILTERED
Assembly k-mer database: $ASSEMBLY_MERYL  
Assembly filtered database: $ASSEMBLY_FILTERED
Reads difference database: $READS_DIFF_ASSEMBLY
Assembly difference database: $ASSEMBLY_DIFF_READS
Intersection database: $INTERSECTION_MERYL

Note: K-mer databases are cached in '$KMER_DB_DIR' for faster subsequent runs.
Analysis uses filtered databases to reduce noise from low-frequency k-mers.
To force regeneration, delete the relevant .meryl directories.

Interpretation:
--------------
Completeness:
The k-mer completeness metric represents the fraction of k-mers from the raw reads
that are also present in the assembly. Higher values indicate
better assembly completeness.

Correctness:
The error rate is calculated based on the k-mers present in the assembly but not in the reads.
The QV score is derived from the error rate, providing a logarithmic scale of quality.
A higher QV indicates better assembly quality.



NOTE: QV score and kmer-completeness metrics ported from merqury. Please refer to the original documentation for more details.
EOF

echo "Results saved to: $RESULTS_FILE"
echo ""
echo "K-mer databases cached in: $KMER_DB_DIR"
echo "Analysis completed using filtered k-mer databases to reduce noise."
echo "Subsequent runs with the same inputs and k-mer size will be faster!"
echo "Analysis complete!"
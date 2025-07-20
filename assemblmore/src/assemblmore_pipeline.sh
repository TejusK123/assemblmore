#!/bin/bash

# assemblmore_pipeline.sh - Main pipeline for genome assembly improvement
# This script orchestrates the complete workflow to improve genome assemblies

set -euo pipefail

# Default values
MAP_PRESET_1="asm20"
MAP_PRESET_2="map-ont"
MAX_ALIGNMENTS=5
OUTPUT_DIR="assemblmore_output"
KEEP_INTERMEDIATE=false
VERBOSE=false
EXPECTED_TELOMERE_LENGTH=8000
LENGTH_THRESHOLD=0
PHRED_THRESHOLD=20
SKIP_BAM=false

usage() {
    cat << EOF
Usage: $0 <reference_genome.fasta> <assembly.fasta> <nanopore_reads.fasta> [OPTIONS]

REQUIRED ARGUMENTS:
  reference_genome.fasta  Path to reference genome (FASTA format)
  assembly.fasta          Path to initial assembly (FASTA format)
  nanopore_reads.fasta    Path to raw sequencing reads (FASTQ/FASTA format)

OPTIONAL ARGUMENTS:
  --expected_telomere_length N    Expected telomere length (default: 8000)
  --phred_threshold N             Phred quality threshold (default: 20)
  --length_threshold N            Minimum contig length threshold (default: 0)
  --output_dir DIR                Output directory (default: assemblmore_output)
  --skip_bam                      Skip BAM file generation (only generate PAF files)
  --help                          Show this help message

DESCRIPTION:
  AssemblMore pipeline for genome assembly improvement through reference-guided
  mapping, gap filling, read alignment, contig spanning, and telomere extension
  with comparative analysis.

PIPELINE STEPS:
  1. Map assembly contigs to reference genome
  2. Create initial refined assembly based on reference mapping
  3. Map reads to refined assembly
  4. Generate final assembly with contig spanning and gap filling
  5. Comparative assembly statistics (original vs improved)

OUTPUT:
  The pipeline generates a comprehensive set of outputs including:
  - Improved assembly: assemblmore_final_assembly.fasta
  - Comparative NX plots and statistics showing improvement
  - Assembly quality metrics and comparison tables
  - Length distribution comparisons
  - Detailed logs and intermediate files

DEPENDENCIES:
  Required: minimap2, samtools, paftools.js, python (with pandas, biopython, etc.)
  Optional: R (with ggplot2, dplyr, readr, scales) for enhanced comparative statistics

EXAMPLES:
  # Basic usage
  $0 reference.fasta assembly.fasta reads.fastq

  # With custom parameters
  $0 reference.fasta assembly.fasta reads.fastq --expected_telomere_length 10000 --phred_threshold 25

  # Skip BAM generation for faster processing
  $0 reference.fasta assembly.fasta reads.fastq --skip_bam

  # Specify output directory
  $0 reference.fasta assembly.fasta reads.fastq --output_dir my_results

EOF
}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Function to log verbose messages
log_verbose() {
    if [ "$VERBOSE" = true ]; then
        log "VERBOSE: $*"
    fi
}

# Function to check if required tools are installed
check_dependencies() {
    local tools=("minimap2" "samtools" "paftools.js" "python" "seqkit")
    local missing=()
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing+=("$tool")
        fi
    done
    
    if [ ${#missing[@]} -ne 0 ]; then
        log "ERROR: Missing required tools: ${missing[*]}"
        log "Please install the missing tools and try again."
        exit 1
    fi
    
    # Check Python dependencies
    local python_deps=("numpy" "pandas" "Bio" "networkx" "more_itertools" "click")
    local missing_python=()
    
    for dep in "${python_deps[@]}"; do
        if ! python -c "import $dep" &> /dev/null; then
            missing_python+=("$dep")
        fi
    done
    
    if [ ${#missing_python[@]} -ne 0 ]; then
        log "ERROR: Missing required Python packages: ${missing_python[*]}"
        log "Install with: pip install numpy pandas biopython networkx more_itertools click"
        exit 1
    fi
    
    # Check R and R packages (optional, but recommended for detailed stats)
    if command -v Rscript &> /dev/null; then
        local r_deps=("ggplot2" "dplyr" "readr" "scales")
        local missing_r=()
        
        for dep in "${r_deps[@]}"; do
            if ! Rscript -e "library($dep)" &> /dev/null; then
                missing_r+=("$dep")
            fi
        done
        
        if [ ${#missing_r[@]} -ne 0 ]; then
            log "WARNING: Missing R packages for detailed assembly statistics: ${missing_r[*]}"
            log "Install with: R -e \"install.packages(c('ggplot2', 'dplyr', 'readr', 'scales'))\""
            log "Will use basic statistics instead"
        else
            log "R and required packages found - will generate detailed assembly statistics"
        fi
    else
        log "WARNING: R not found - will use basic assembly statistics only"
        log "Install R and packages for comprehensive analysis"
    fi
    
    log "All required dependencies found"
}

# Function to validate input files
validate_inputs() {
    local ref_genome="$1"
    local assembly="$2"
    local reads="$3"
    
    if [ ! -f "$ref_genome" ]; then
        log "ERROR: Reference genome file '$ref_genome' does not exist"
        exit 1
    fi
    
    if [ ! -f "$assembly" ]; then
        log "ERROR: Assembly file '$assembly' does not exist"
        exit 1
    fi
    
    if [ ! -f "$reads" ]; then
        log "ERROR: Nanopore reads file '$reads' does not exist"
        exit 1
    fi
    
    log "All input files validated"
}

# Function to get absolute path
get_abs_path() {
    local input_path="$1"
    
    # If already absolute path, return as-is
    if [[ "$input_path" = /* ]]; then
        echo "$input_path"
        return 0
    fi
    
    # For relative paths, use realpath if available, otherwise manual resolution
    if command -v realpath &> /dev/null; then
        realpath "$input_path" 2>/dev/null || {
            # Fallback to manual resolution if realpath fails
            echo "$(pwd)/$input_path"
        }
    else
        # Manual resolution for systems without realpath
        local dir_part=$(dirname "$input_path")
        local file_part=$(basename "$input_path")
        
        # Try to cd to directory and get absolute path
        if [ -d "$dir_part" ]; then
            echo "$(cd "$dir_part" && pwd)/$file_part"
        else
            # If directory doesn't exist, return the path as constructed from pwd
            echo "$(pwd)/$input_path"
        fi
    fi
}

# Function to cleanup intermediate files
cleanup() {
    if [ "$KEEP_INTERMEDIATE" = false ]; then
        log "Cleaning up intermediate files..."
        # Remove intermediate files but keep final outputs
        find "$OUTPUT_DIR" -name "*.sam" -delete 2>/dev/null || true
        find "$OUTPUT_DIR" -name "*_mapped_to_*.sorted.paf" -not -name "final_*" -delete 2>/dev/null || true
        # Only remove BAM files if they were generated (i.e., SKIP_BAM is false)
        if [ "$SKIP_BAM" = false ]; then
            find "$OUTPUT_DIR" -name "*_mapped_to_*.sorted.bam*" -delete 2>/dev/null || true
        fi
    else
        log "Keeping all intermediate files as requested"
    fi
}

# Function to run a step with error handling
run_step() {
    local step_name="$1"
    local step_cmd="$2"
    
    log "Starting $step_name..."
    log_verbose "Command: $step_cmd"
    
    if eval "$step_cmd"; then
        log "✓ $step_name completed successfully"
    else
        log "✗ $step_name failed"
        exit 1
    fi
}

# Parse command line arguments
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir|--output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--preset1)
            MAP_PRESET_1="$2"
            shift 2
            ;;
        -P|--preset2)
            MAP_PRESET_2="$2"
            shift 2
            ;;
        -N|--max-alignments)
            MAX_ALIGNMENTS="$2"
            shift 2
            ;;
        -t|--telomere-length|--expected_telomere_length)
            EXPECTED_TELOMERE_LENGTH="$2"
            shift 2
            ;;
        -l|--length-threshold|--length_threshold)
            LENGTH_THRESHOLD="$2"
            shift 2
            ;;
        -q|--phred-threshold|--phred_threshold)
            PHRED_THRESHOLD="$2"
            shift 2
            ;;
        -k|--keep-intermediate)
            KEEP_INTERMEDIATE=true
            shift
            ;;
        --skip-bam|--skip_bam)
            SKIP_BAM=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            log "ERROR: Unknown option $1"
            usage
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Check for required positional arguments
if [ ${#POSITIONAL_ARGS[@]} -ne 3 ]; then
    log "ERROR: Exactly 3 positional arguments required <reference_genome.fasta> <assembly.fasta> <nanopore_reads.fasta>"
    log "Provided: ${POSITIONAL_ARGS[*]}"
    usage
    exit 1
fi

REFERENCE_GENOME="${POSITIONAL_ARGS[0]}"
ASSEMBLY="${POSITIONAL_ARGS[1]}"
NANOPORE_READS="${POSITIONAL_ARGS[2]}"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Main pipeline execution
main() {
    log "Starting assemblmore pipeline"
    log "Reference genome: $REFERENCE_GENOME"
    log "Assembly: $ASSEMBLY"
    log "Nanopore reads: $NANOPORE_READS"
    log "Output directory: $OUTPUT_DIR"
    
    # Check dependencies
    check_dependencies
    
    # Get absolute paths for input files (before changing directories)
    REF_ABS=$(get_abs_path "$REFERENCE_GENOME")
    ASM_ABS=$(get_abs_path "$ASSEMBLY")
    READS_ABS=$(get_abs_path "$NANOPORE_READS")
    
    # Validate inputs using absolute paths
    validate_inputs "$REF_ABS" "$ASM_ABS" "$READS_ABS"
    
    # Create output directory and change to it
    mkdir -p "$OUTPUT_DIR"
    cd "$OUTPUT_DIR"
    
    # Generate base names for file tracking
    REF_BASE=$(basename "$REF_ABS" | sed 's/\.\(fasta\|fa\|fna\)$//')
    ASM_BASE=$(basename "$ASM_ABS" | sed 's/\.\(fasta\|fa\|fna\)$//')
    READS_BASE=$(basename "$READS_ABS" | sed 's/\.\(fastq\|fq\|fasta\|fa\)$//')


    REF_BASE="${REF_BASE%.fasta}"
    REF_BASE="${REF_BASE%.fa}"
    REF_BASE="${REF_BASE%.fna}"
    ASM_BASE="${ASM_BASE%.fasta}"
    ASM_BASE="${ASM_BASE%.fa}"
    ASM_BASE="${ASM_BASE%.fna}"
    READS_BASE="${READS_BASE%.fastq}"
    READS_BASE="${READS_BASE%.fq}"
    READS_BASE="${READS_BASE%.fna}"
    READS_BASE="${READS_BASE%.fasta}"
    READS_BASE="${READS_BASE%.fa}"
    
    log "$REF_BASE"
    log "$ASM_BASE"
    log "$READS_BASE"

    log_verbose "Reference base name: $REF_BASE"
    log_verbose "Assembly base name: $ASM_BASE"
    log_verbose "Reads base name: $READS_BASE"
    log_verbose "Skip BAM generation: $SKIP_BAM"
    #log_verbose "PAF file: $STEP1_PAF"
    #log_verbose "Actual reference base (from PAF): $ACTUAL_REF_BASE"
    log_verbose "Expected telomere length: $EXPECTED_TELOMERE_LENGTH"
    log_verbose "Length threshold: $LENGTH_THRESHOLD"
    log_verbose "Phred threshold: $PHRED_THRESHOLD"
    
    # Step 1: Map assembly contigs to reference genome



    STEP1_PAF="${ASM_BASE}_mapped_to_${REF_BASE}.sorted.paf"
    BAM_FLAG=""
    if [ "$SKIP_BAM" = true ]; then
        BAM_FLAG="--no-bam"
    fi

    if [ ! -f "$STEP1_PAF" ]; then
        run_step "Step 1: Mapping assembly contigs to reference" \
            "\"$SCRIPT_DIR/fill_gaps.sh\" $BAM_FLAG \"$REF_ABS\" \"$ASM_ABS\" \"$MAP_PRESET_1\""
    else
        log "Skipping Step 1: PAF file already exists: $STEP1_PAF"
    fi

    if [ ! -f "$STEP1_PAF" ]; then
        log "ERROR: Step 1 did not produce expected output file: $STEP1_PAF"
        exit 1
    fi
    
    # Step 2: Determine contig placements and create initial refined assembly
    # Note: The actual output names are determined by contig_placements.py based on PAF filename parsing
    # We need to figure out what the actual output names will be
    PAF_BASE=$(basename "$STEP1_PAF" .sorted.paf)
    if [[ "$PAF_BASE" == *"_mapped_to_"* ]]; then
        ACTUAL_REF_BASE="${PAF_BASE##*_mapped_to_}"
    else
        ACTUAL_REF_BASE="$REF_BASE"
    fi
    
    ACTUAL_REF_BASE="${ACTUAL_REF_BASE%.fasta}"
    ACTUAL_REF_BASE="${ACTUAL_REF_BASE%.fa}"
    ACTUAL_REF_BASE="${ACTUAL_REF_BASE%.fna}"

    FILTERED_TSV="filtered_by_${ACTUAL_REF_BASE}_contigs.tsv"
    INITIAL_ASSEMBLY="ordered_and_oriented_to_${ACTUAL_REF_BASE}_assembly.fasta"
    run_step "Step 2: Creating initial refined assembly" \
        "\"$SCRIPT_DIR/initial_assembly.sh\" \"$STEP1_PAF\" \"$READS_ABS\" \"$ASM_ABS\""
    
    if [ ! -f "$FILTERED_TSV" ] || [ ! -f "$INITIAL_ASSEMBLY" ]; then
        log "ERROR: Step 2 did not produce expected output files"
        log "Expected: $FILTERED_TSV and $INITIAL_ASSEMBLY"
        exit 1
    fi
    
    # Step 3: Map reads to the partially refined assembly
    INITIAL_BASE="ordered_and_oriented_to_${ACTUAL_REF_BASE}_assembly"
    STEP3_PAF="${READS_BASE}_mapped_to_${INITIAL_BASE}.sorted.paf"
    
    if [ ! -f "$STEP3_PAF" ]; then
    run_step "Step 3: Mapping reads to refined assembly" \
        "\"$SCRIPT_DIR/fill_gaps.sh\" $BAM_FLAG \"$INITIAL_ASSEMBLY\" \"$READS_ABS\" \"$MAP_PRESET_2\" \"$MAX_ALIGNMENTS\""
    else
        log "Skipping Step 3: PAF file already exists: $STEP3_PAF"
    fi
    
    if [ ! -f "$STEP3_PAF" ]; then
        log "ERROR: Step 3 did not produce expected output file: $STEP3_PAF"
        exit 1
    fi
    
    # Step 4: Generate final assembly
    run_step "Step 4: Generating final assembly" \
        "\"$SCRIPT_DIR/span_contigs.sh\" \"$STEP3_PAF\" \"$FILTERED_TSV\" \"$INITIAL_ASSEMBLY\" \"$READS_ABS\" --expected_telomere_length \"$EXPECTED_TELOMERE_LENGTH\" --length_threshold \"$LENGTH_THRESHOLD\" --phred_threshold \"$PHRED_THRESHOLD\""
    
    # Rename final output to something more descriptive
    FINAL_OUTPUT="assemblmore_final_assembly.fasta"
    if [ -f "final_assembly.fasta" ]; then
        mv "final_assembly.fasta" "$FINAL_OUTPUT"
        log "✓ Final assembly saved as: $FINAL_OUTPUT"
    else
        # Look for other potential output files from span_contigs.py
        for possible_output in *.fasta; do
            if [[ "$possible_output" != "$INITIAL_ASSEMBLY" && "$possible_output" != "assemblmore_final_assembly.fasta" ]]; then
                mv "$possible_output" "$FINAL_OUTPUT"
                log "✓ Final assembly saved as: $FINAL_OUTPUT"
                break
            fi
        done
    fi
    
    # Cleanup intermediate files if requested
    cleanup
    
    log "✓ Pipeline completed successfully!"
    log "Output files are located in: $(pwd)"
    log "Key outputs:"
    log "  - Contig placements: $FILTERED_TSV"
    log "  - Initial refined assembly: $INITIAL_ASSEMBLY"
    log "  - Final assembly: $FINAL_OUTPUT"
    
    # Generate comprehensive assembly statistics with comparison
    if command -v Rscript &> /dev/null; then
        log ""
        log "Generating comprehensive assembly statistics with comparison..."
        
        # Run custom R script for comparative assembly analysis
        run_step "Generating comparative assembly statistics" \
            "Rscript \"$SCRIPT_DIR/assembly_stats.R\" \"$ASM_ABS:Original\" \"$(pwd)/$FINAL_OUTPUT:Improved\" \"$(pwd)/stats\""

        log "Comparative assembly statistics saved in: $(pwd)/stats/"
        log "  - Individual NX plots: stats/*_nx_plot.png"
        log "  - Comparative NX plot: stats/comparative_nx_plot.png"
        log "  - Assembly comparison table: stats/assembly_comparison.csv"
        log "  - Length distribution comparison: stats/length_distribution_comparison.png"
        log "  - Comprehensive summary: stats/assembly_analysis_summary.txt"
    else
        log "Rscript not found - skipping detailed assembly statistics"
        log "Install R and required packages (ggplot2, dplyr, readr, scales) for comprehensive analysis"
        
        # Fallback to seqkit if available
        if command -v seqkit &> /dev/null; then
            log ""
            log "Basic assembly statistics (seqkit):"
            log "Original assembly:"
            seqkit stats "$ASM_ABS" | tail -n +2
            log "Final assembly:"
            seqkit stats "$FINAL_OUTPUT" | tail -n +2
        fi
    fi
}

# Set up trap for cleanup on exit
trap 'cleanup' EXIT

# Run main pipeline
main "$@"

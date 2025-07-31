#!/bin/bash

# AssemblMore Example: C. briggsae AF16 Strain
# 
# This script demonstrates AssemblMore pipeline using:
# - Organism: Caenorhabditis briggsae AF16 strain  
# - Sequencing: MinION r10.4 Nanopore (Rog Lab)
# - Reference: C. briggsae QX1410 (PRJNA784955, WBPS19)
# - Assembly: Flye-generated assembly from long reads

set -euo pipefail

EXAMPLE_DIR="assemblmore_example_briggsae_AF16"
DOWNLOAD_URL="https://drive.google.com/drive/folders/1b03BW8kaFKKbbPbfvOPKR6ZUZ5fskX6a?usp=drive_link"

echo "=== AssemblMore C. briggsae AF16 Example ==="
echo "This example demonstrates genome assembly improvement using:"
echo "  - Reference: C. briggsae QX1410 genome"
echo "  - Assembly: Flye assembly of AF16 strain"
echo "  - Reads: MinION r10.4 Nanopore sequencing data"
echo ""

# Check if gdown is installed
if ! command -v gdown &> /dev/null; then
    echo "ERROR: gdown is required to download example data"
    echo "Install with: pip install gdown"
    exit 1
fi

# Check if assemblmore is available
if ! command -v assemblmore &> /dev/null; then
    echo "ERROR: assemblmore command not found"
    echo "Make sure AssemblMore src directory is in your PATH"
    echo "See installation instructions in README.md"
    exit 1
fi

# Download and prepare data
if [ -d "$EXAMPLE_DIR" ]; then
    echo "Example data directory already exists: $EXAMPLE_DIR"
else
    echo "Downloading example data from Google Drive..."
    echo "URL: $DOWNLOAD_URL"
    gdown --folder "$DOWNLOAD_URL"
    
    if [ -d "$EXAMPLE_DIR" ]; then
        echo "Decompressing files..."
        cd "$EXAMPLE_DIR" && gunzip -d *.fasta.gz *.fa.gz 2>/dev/null || true
        cd ..
        echo "✓ Data downloaded and prepared"
    else
        echo "ERROR: Failed to download example data"
        echo "Please check your internet connection and try again"
        exit 1
    fi
fi

# Verify required files exist
REFERENCE="$EXAMPLE_DIR/caenorhabditis_briggsae.QX1410_PRJNA784955.WBPS19.genomic.fa"
ASSEMBLY="$EXAMPLE_DIR/assembly_AF16_flye.fasta"
READS="$EXAMPLE_DIR/C_briggsae_AF16.fasta"

for file in "$REFERENCE" "$ASSEMBLY" "$READS"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file not found: $file"
        echo "Please check the download completed successfully"
        exit 1
    fi
done

echo ""
echo "=== Running AssemblMore Pipeline ==="
echo "Reference: $REFERENCE"
echo "Assembly:  $ASSEMBLY" 
echo "Reads:     $READS"
echo ""
echo "Starting pipeline... (this may take 30-60 minutes)"

# Run AssemblMore pipeline
assemblmore "$REFERENCE" "$ASSEMBLY" "$READS" --output_dir "assemblmore_output_example"

echo ""
echo "=== Pipeline Completed Successfully! ==="
echo ""
echo "Output files located in: assemblmore_output_example/"
echo ""
echo "Key results:"
echo "  - Final assembly: assemblmore_output_example/assemblmore_final_assembly.fasta"
echo "  - Coverage data:  assemblmore_output_example/final_assembly_coverage.txt"
echo "  - Statistics:     assemblmore_output_example/stats/"
echo ""
echo "Check the following for assembly improvements:"
echo "  - assemblmore_output_example/stats/comparative_nx_plot.png"
echo "  - assemblmore_output_example/stats/assembly_comparison.csv"
echo "  - assemblmore_output_example/stats/assembly_analysis_summary.txt"
echo ""
echo "For more information, see: assemblmore/examples/README.md"
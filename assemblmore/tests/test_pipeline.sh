#!/bin/bash

# test_pipeline.sh - Test script for AssemblMore pipeline

set -euo pipefail

echo "Testing AssemblMore pipeline..."

# Create test directory
TEST_DIR="test_assemblmore"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

# Function to create a simple test FASTA file
create_test_fasta() {
    local filename="$1"
    local num_seqs="$2"
    
    cat > "$filename" << EOF
>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
EOF
}

echo "Creating test files..."
create_test_fasta "reference.fasta" 2
create_test_fasta "assembly.fasta" 2
create_test_fasta "reads.fasta" 2

echo "Testing dependency check..."
if ../src/assemblmore_pipeline.sh --help > /dev/null 2>&1; then
    echo "✓ Pipeline script is executable and shows help"
else
    echo "✗ Pipeline script failed"
    exit 1
fi

echo "Testing with invalid arguments..."
if ../src/assemblmore_pipeline.sh 2>/dev/null; then
    echo "✗ Pipeline should fail with no arguments"
    exit 1
else
    echo "✓ Pipeline correctly rejects invalid arguments"
fi

echo "Test files created in: $(pwd)"
echo "To test the full pipeline (requires real data):"
echo "  ../src/assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fasta"

cd ..
echo "Test completed successfully!"

#!/bin/bash

# test_arguments.sh - Test script to verify parameter passing

echo "Testing assemblmore pipeline argument parsing..."

echo ""
echo "=== Testing help output ==="
./assemblmore_pipeline.sh --help | head -20

echo ""
echo "=== Testing with custom parameters ==="
echo "Command that would be executed:"
echo "./assemblmore_pipeline.sh -t 10000 -l 5000 -q 30 -v ref.fasta asm.fasta reads.fasta"

echo ""
echo "=== Testing argument validation ==="
if ./assemblmore_pipeline.sh -t 10000 2>&1 | grep -q "ERROR: Exactly 3 positional arguments required"; then
    echo "✓ Argument validation works correctly"
else
    echo "✗ Argument validation failed"
fi

echo ""
echo "=== Testing unknown option handling ==="
if ./assemblmore_pipeline.sh --unknown-option 2>&1 | grep -q "ERROR: Unknown option"; then
    echo "✓ Unknown option handling works correctly"
else
    echo "✗ Unknown option handling failed"
fi

echo ""
echo "All argument parsing tests completed!"

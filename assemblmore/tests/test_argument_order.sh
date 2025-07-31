#!/bin/bash

# test_argument_order.sh - Test argument parsing with different orders

echo "Testing argument order flexibility..."

echo ""
echo "=== Test 1: Options before positional arguments ==="
echo "Command: ../src/assemblmore_pipeline.sh -t 10000 -v ref.fasta asm.fasta reads.fasta"
if ../src/assemblmore_pipeline.sh -t 10000 -v ref.fasta asm.fasta reads.fasta 2>&1 | grep -q "Reference genome: ref.fasta"; then
    echo "✓ Options before positional works"
else
    echo "✗ Options before positional failed"
fi

echo ""
echo "=== Test 2: Options after positional arguments ==="
echo "Command: ../src/assemblmore_pipeline.sh ref.fasta asm.fasta reads.fasta -t 10000 -v"
if ../src/assemblmore_pipeline.sh ref.fasta asm.fasta reads.fasta -t 10000 -v 2>&1 | grep -q "Reference genome: ref.fasta"; then
    echo "✓ Options after positional works"
else
    echo "✗ Options after positional failed"
fi

echo ""
echo "=== Test 3: Mixed order ==="
echo "Command: ../src/assemblmore_pipeline.sh -t 10000 ref.fasta -v asm.fasta reads.fasta -q 30"
if ../src/assemblmore_pipeline.sh -t 10000 ref.fasta -v asm.fasta reads.fasta -q 30 2>&1 | grep -q "Reference genome: ref.fasta"; then
    echo "✓ Mixed order works"
else
    echo "✗ Mixed order failed"
fi

echo ""
echo "=== Test 4: Wrong number of positional arguments ==="
echo "Command: ../src/assemblmore_pipeline.sh -t 10000 ref.fasta asm.fasta"
if ../src/assemblmore_pipeline.sh -t 10000 ref.fasta asm.fasta 2>&1 | grep -q "Exactly 3 positional arguments required"; then
    echo "✓ Correctly rejects wrong number of arguments"
else
    echo "✗ Failed to reject wrong number of arguments"
fi

echo ""
echo "All argument order tests completed!"

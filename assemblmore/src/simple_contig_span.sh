#!/bin/bash


conda activate Nanopore
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$SCRIPT_DIR/simple_contig_span_v2.py" "$@"

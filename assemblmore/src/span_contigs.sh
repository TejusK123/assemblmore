#!/bin/bash



set -euo pipefail
# Check if Python is installed

if ! command -v python &> /dev/null; then
    echo "Error: python is not installed or not in PATH." >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python "$SCRIPT_DIR/span_contigs.py" "$@"
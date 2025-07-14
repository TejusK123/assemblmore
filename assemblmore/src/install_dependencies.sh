#!/bin/bash

# install_dependencies.sh - Install dependencies for AssemblMore

set -euo pipefail

echo "Installing dependencies for AssemblMore..."

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for package managers
if command_exists conda; then
    echo "Using conda to install dependencies..."
    conda install -c bioconda minimap2 samtools seqkit merqury -y
    pip install numpy pandas biopython networkx more_itertools click
elif command_exists brew; then
    echo "Using homebrew to install dependencies..."
    brew install minimap2 samtools seqkit
    pip install numpy pandas biopython networkx more_itertools click
else
    echo "Please install the following tools manually:"
    echo "  - minimap2"
    echo "  - samtools" 
    echo "  - seqkit"
    echo ""
    echo "And install Python packages with:"
    echo "  pip install numpy pandas biopython networkx more_itertools click"
    exit 1
fi

echo "Dependencies installed successfully!"
echo "You can now run AssemblMore with: ./assemblmore_pipeline.sh"

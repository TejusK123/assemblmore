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
    # Core bioinformatics tools
    conda install -c bioconda -c conda-forge minimap2 samtools seqkit merqury k8 meryl
    # Python packages
    pip install numpy pandas biopython networkx more_itertools click
    # Optional R packages for enhanced statistics
    echo "Installing optional R packages for enhanced statistics..."
    if command_exists R; then
        R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
        R -e "install.packages(c('ggplot2', 'dplyr', 'readr', 'scales', 'dbscan'))"
    else
        echo "R not found - skipping R package installation"
    fi
elif command_exists brew; then
    echo "Using homebrew to install dependencies..."
    # Core bioinformatics tools
    brew install minimap2 samtools seqkit
    # Python packages
    pip install numpy pandas biopython networkx more_itertools click
    # Try to install additional tools that may be available via brew
    if brew list meryl &>/dev/null || brew install meryl; then
        echo "meryl installed via homebrew"
    else
        echo "Warning: meryl not available via homebrew - k-mer completeness analysis will be disabled"
    fi
    # Optional R packages
    echo "Installing optional R packages for enhanced statistics..."
    if command_exists R; then
        R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
        R -e "install.packages(c('ggplot2', 'dplyr', 'readr', 'scales', 'dbscan'), dependencies = TRUE)"
    else
        echo "R not found - install R for enhanced statistical analysis and visualization"
    fi
else
    echo "Please install the following tools manually:"
    echo ""
    echo "REQUIRED TOOLS:"
    echo "  - minimap2 (sequence alignment)"
    echo "  - samtools (SAM/BAM manipulation)" 
    echo "  - seqkit (FASTA/FASTQ manipulation)"
    echo "  - k8 (JavaScript engine for paftools.js - usually comes with minimap2)"
    echo ""
    echo "OPTIONAL TOOLS:"
    echo "  - meryl (for k-mer completeness analysis)"
    echo "  - merqury (for assembly quality assessment)"
    echo "  - R (for enhanced statistics and visualization)"
    echo ""
    echo "PYTHON PACKAGES (required):"
    echo "  pip install numpy pandas biopython networkx more_itertools click"
    echo ""
    echo "R PACKAGES (optional, for enhanced analysis):"
    echo "  R -e \"install.packages(c('ggplot2', 'dplyr', 'readr', 'scales', 'dbscan'))\""
    echo ""
    echo "Note: paftools.js is typically included with minimap2 installation"
    echo "      If not available, ensure k8 JavaScript engine is installed"
    exit 1
fi

echo ""
echo "=== INSTALLATION COMPLETE ==="
echo ""
echo "Core dependencies installed successfully!"
echo ""
echo "INSTALLED TOOLS:"
echo "  ✓ minimap2 - sequence alignment"
echo "  ✓ samtools - SAM/BAM manipulation"
echo "  ✓ seqkit - FASTA/FASTQ manipulation"
echo "  ✓ Python packages (numpy, pandas, biopython, networkx, more_itertools, click)"

# Check for optional tools
if command_exists meryl; then
    echo "  ✓ meryl - k-mer analysis (optional)"
else
    echo "  ⚠ meryl - not installed (k-mer completeness analysis will be disabled)"
fi

if command_exists merqury; then
    echo "  ✓ merqury - assembly quality assessment (optional)"
else
    echo "  ⚠ merqury - not installed (some quality metrics may be unavailable)"
fi

if command_exists R; then
    echo "  ✓ R - statistical analysis and visualization (optional)"
else
    echo "  ⚠ R - not installed (enhanced statistics and plots will be unavailable)"
fi

if command_exists k8; then
    echo "  ✓ k8 - JavaScript engine for paftools.js"
else
    echo "  ⚠ k8 - not found (paftools.js may not work properly)"
fi

echo ""
echo "You can now run AssemblMore with:"
echo "  ./assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fastq"
echo ""
echo "For help and usage examples:"
echo "  ./assemblmore_pipeline.sh --help"

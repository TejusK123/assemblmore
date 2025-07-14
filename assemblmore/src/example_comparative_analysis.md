# AssemblMore Comparative Assembly Analysis

## Overview

The AssemblMore pipeline now includes comprehensive comparative analysis capabilities that allow you to directly compare multiple assemblies in a single analysis.

## Key Features

### 🔍 **Comparative Statistics**
- Side-by-side comparison of key metrics (N50, auN, contig count, etc.)
- Automated generation of comparison tables in CSV format
- Statistical summaries for multiple assemblies

### 📊 **Visual Comparisons**
- **Comparative NX plots** - overlay NX curves for direct visual comparison
- **Length distribution comparisons** - histograms showing contig size distributions
- **Individual assembly plots** - detailed plots for each assembly

### 📈 **Advanced Metrics**
- **auN calculation** - area under the NX curve for single-value quality assessment
- **NX values** - N10, N50, N90, and full NX spectrum
- **Basic statistics** - total length, contig count, mean/median lengths

## Usage Examples

### Single Assembly Analysis
```bash
# Basic analysis of one assembly
Rscript assembly_stats.R assembly.fasta:MyAssembly output_dir
```

### Comparative Analysis (Pipeline Integration)
```bash
# The pipeline automatically compares original vs improved assembly
./assemblmore_pipeline.sh reads.fastq assembly.fasta --output_dir results
```

### Manual Comparative Analysis
```bash
# Compare multiple assemblies manually
Rscript assembly_stats.R \
    original.fasta:Original \
    improved.fasta:Improved \
    polished.fasta:Polished \
    results_dir
```

### Using Default Names
```bash
# Assembly names derived from filenames
Rscript assembly_stats.R assembly1.fasta assembly2.fasta results_dir
```

## Output Files

### Individual Assembly Files
- `{assembly_name}_nx_plot.png` - Individual NX curve plot
- `{assembly_name}_nx_values.csv` - NX values data
- Individual statistics in summary report

### Comparative Files (Multiple Assemblies)
- `comparative_nx_plot.png` - Overlay NX curves for all assemblies
- `assembly_comparison.csv` - Side-by-side comparison table
- `length_distribution_comparison.png` - Contig length distributions
- `assembly_analysis_summary.txt` - Comprehensive summary report

## Integration with Pipeline

The AssemblMore pipeline automatically generates comparative analysis between:
- **Original assembly** - your input assembly
- **Improved assembly** - the pipeline's improved output

This allows you to immediately see the quantitative improvements achieved by the pipeline, including:
- Changes in N50 values
- Reduction in contig count (due to gap filling and spanning)
- Improvements in contiguity metrics
- Visual comparison of assembly quality

## Dependencies

### Required (Core Functionality)
- R with packages: ggplot2, dplyr, readr, scales

### Installation
```bash
# Install R packages
R -e "install.packages(c('ggplot2', 'dplyr', 'readr', 'scales'))"
```

## Benefits

1. **Immediate Quality Assessment** - See improvements at a glance
2. **Publication-Ready Figures** - High-quality NX plots and comparisons
3. **Quantitative Metrics** - Precise measurements of assembly quality
4. **Streamlined Workflow** - Integrated into the pipeline for automatic comparison
5. **Flexible Usage** - Can compare any number of assemblies manually

This enhancement makes AssemblMore not just an assembly improvement tool, but also a comprehensive assembly analysis and comparison platform.

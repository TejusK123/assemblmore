# AssemblMore Examples

This directory contains example datasets and scripts to demonstrate the AssemblMore pipeline functionality.

## Available Examples

### C. briggsae AF16 Strain Example

**Dataset**: Caenorhabditis briggsae AF16 strain  
**Source**: Rog Lab  
**Sequencing**: MinION r10.4 Nanopore sequencing  
**Species**: Caenorhabditis briggsae (nematode)  

This example demonstrates AssemblMore's ability to improve genome assemblies using:
- **Reference genome**: C. briggsae QX1410 reference (PRJNA784955, WBPS19)
- **Initial assembly**: Flye assembly of AF16 strain long reads
- **Raw reads**: MinION r10.4 Nanopore sequencing data

## Quick Start

### Prerequisites

1. **Install AssemblMore dependencies**:
   ```bash
   cd ../src
   ./install_dependencies.sh
   ```

2. **Install gdown** (for downloading example data):
   ```bash
   pip install gdown
   ```

### Running the Example

```bash
# Navigate to examples directory
cd assemblmore/examples

# Make the script executable
chmod +x run_example.sh

# Run the example
./run_example.sh
```

### What the Example Does

1. **Downloads data** (if not already present):
   - Reference genome: `caenorhabditis_briggsae.QX1410_PRJNA784955.WBPS19.genomic.fa`
   - Initial assembly: `assembly_AF16_flye.fasta`
   - Raw reads: `C_briggsae_AF16.fasta`

2. **Runs AssemblMore pipeline**:
   - Maps assembly contigs to reference genome
   - Performs contig placement and orientation
   - Fills gaps using spanning reads
   - Generates improved assembly with telomere extension
   - Produces comparative statistics and visualizations

### Expected Output

The pipeline will create an `assemblmore_output` directory containing:

#### Primary Results
- `assemblmore_final_assembly.fasta` - Improved assembly
- `final_assembly_coverage.txt` - Coverage analysis
- `assembly_comparison.csv` - Comparative statistics

#### Analysis and Visualization
- `stats/comparative_nx_plot.png` - NX curve comparison
- `stats/length_distribution_comparison.png` - Contig length distributions
- `stats/assembly_comparison.csv` - Detailed metrics comparison
- `stats/assembly_analysis_summary.txt` - Comprehensive summary

#### Intermediate Files
- `filtered_by_*_contigs.tsv` - Contig placement results
- `ordered_and_oriented_*.fasta` - Initial refined assembly
- Various alignment files (PAF/BAM format)

### Performance Metrics

The C. briggsae AF16 example typically shows:
- **Improved contiguity**: Higher N50, reduced contig count
- **Better orientation**: Corrected contig placement relative to reference
- **Gap filling**: Reduced number of gaps and ambiguous bases
- **Enhanced coverage**: More uniform coverage distribution

### Customizing Parameters

You can modify the example script to test different parameters:

```bash
# For faster processing (skip BAM generation)
assemblmore reference.fa assembly.fa reads.fa --skip_bam

# Custom telomere length
assemblmore reference.fa assembly.fa reads.fa --expected_telomere_length 12000

# Higher quality threshold
assemblmore reference.fa assembly.fa reads.fa --phred_threshold 25

# Custom output directory
assemblmore reference.fa assembly.fa reads.fa --output_dir my_results
```

## Dataset Details

### C. briggsae AF16 Strain

- **Organism**: *Caenorhabditis briggsae*
- **Strain**: AF16
- **Genome size**: ~104 Mb
- **Chromosomes**: 6 (I, II, III, IV, V, X)
- **Sequencing platform**: Oxford Nanopore MinION
- **Chemistry**: r10.4 flow cells
- **Assembly method**: Flye assembler

### Reference Genome

- **Assembly**: QX1410_PRJNA784955
- **Source**: WormBase ParaSite (WBPS19)
- **Quality**: High-quality reference with chromosome-level assembly
- **Annotation**: Comprehensive gene annotations available

## Troubleshooting

### Common Issues

1. **Download failures**:
   ```bash
   # Manual download if gdown fails
   # Visit: https://drive.google.com/drive/folders/1b03BW8kaFKKbbPbfvOPKR6ZUZ5fskX6a
   ```

2. **Memory issues**:
   ```bash
   # Use --skip_bam for reduced memory usage
   ./run_example.sh --skip_bam
   ```

3. **Missing dependencies**:
   ```bash
   # Check and install missing tools
   cd ../src && ./install_dependencies.sh
   ```

### Getting Help

- **Pipeline help**: `assemblmore --help`
- **GitHub Issues**: [Report problems](https://github.com/TejusK123/assemblmore/issues)
- **Documentation**: See main README.md for detailed usage

## Contributing Examples

To contribute additional examples:

1. **Prepare dataset**: Ensure data is publicly accessible
2. **Create script**: Follow the pattern in `run_example.sh`
3. **Document thoroughly**: Include organism info, sequencing details
4. **Test thoroughly**: Verify the example works on different systems
5. **Submit PR**: Include documentation updates

### Example Script Template

```bash
#!/bin/bash

# [ORGANISM] [STRAIN] Example
# Source: [DATA_SOURCE]
# Sequencing: [PLATFORM] [CHEMISTRY]

EXAMPLE_DIR="assemblmore_example_[organism]_[strain]"
DOWNLOAD_URL="[YOUR_DOWNLOAD_URL]"

if [ -d "$EXAMPLE_DIR" ]; then
    echo "Example data already downloaded."
else
    echo "Downloading example data..."
    # Your download method here
fi

assemblmore $EXAMPLE_DIR/reference.fa \\
           $EXAMPLE_DIR/assembly.fa \\
           $EXAMPLE_DIR/reads.fa
```

## Acknowledgments

- **Rog Lab** for providing the C. briggsae AF16 dataset
- **WormBase ParaSite** for reference genome data
- **Oxford Nanopore Technologies** for sequencing platform

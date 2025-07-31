# AssemblMore
A comprehensive genome assembly improvement pipeline that automatically merges and refines contigs from next-generation assemblers (Flye, Canu, etc.) using reference-guided approaches and ultra-long reads.

![GitHub](https://img.shields.io/github/license/TejusK123/assemblmore)
![GitHub last commit](https://img.shields.io/github/last-commit/TejusK123/assemblmore)

## Overview

AssemblMore enhances genome assemblies through a multi-step pipeline that:
- Maps assembly contigs to reference genomes for optimal placement
- Uses spanning ultra-long reads to merge and extend contigs
- Fills gaps with high-quality sequence data
- Provides comprehensive quality assessment and comparative analysis
- Future iterations will attempt to 'un-collapse' complex genomic regions such as rRNA arrays, telomeres, etc.

## Features

### 🧬 **Assembly Improvement**
- Reference-guided contig placement and orientation
- Ultra-long read spanning for contig merging
- Intelligent gap filling with quality filtering
- Telomere extension and complex region handling

### 📊 **Quality Assessment**
- Comprehensive assembly statistics (N50, auN, contig count)
- Visual comparative analysis with NX plots
- Coverage analysis and quality metrics
- Before/after improvement comparisons

### 🔧 **Robust Pipeline**
- Automated dependency checking and installation
- Flexible parameter configuration
- Comprehensive error handling and logging
- Support for various input formats (FASTA/FASTQ)

### 📈 **Advanced Analytics**
- Issues metric calculation with cluster identification
- Length distribution analysis
- Coverage depth assessment
- R-based statistical visualization

## Quick Start

### Installation

1. **Clone the repository:**
```bash
git clone https://github.com/TejusK123/assemblmore.git
cd assemblmore
```

2. **Install dependencies:**
```bash
cd assemblmore/src
chmod +x install_dependencies.sh
./install_dependencies.sh
```

3. **Add to PATH (optional):**
```bash
# Add to your ~/.bashrc or ~/.zshrc
export PATH="$HOME/path/to/assemblmore/assemblmore/src:$PATH"
```

### Basic Usage

```bash

# Run the pipeline
<path>/assemblmore reference.fasta assembly.fasta reads.fastq
```

## Installation

### Prerequisites

The pipeline requires the following tools:
- **minimap2** - For sequence alignment
- **samtools** - For SAM/BAM file manipulation  
- **paftools.js** - For PAF format conversion
- **seqkit** - For FASTA/FASTQ manipulation
- **python** (3.6+) - For analysis scripts
- **R** (optional) - For enhanced statistics and visualization

### Python Dependencies
```bash
pip install numpy pandas biopython networkx more_itertools click
```

### R Dependencies (Optional)
```r
install.packages(c("ggplot2", "dplyr", "readr", "scales", "dbscan"))
```

### Automated Installation

Use the provided script to install all dependencies:
```bash
cd assemblmore/src
./install_dependencies.sh
```

## Usage

### Command Line Interface

```bash
./assemblmore_pipeline.sh <reference_genome.fasta> <assembly.fasta> <reads.fastq> [OPTIONS]
```

### Required Arguments
- `reference_genome.fasta` - Reference genome in FASTA format
- `assembly.fasta` - Initial assembly to improve
- `reads.fastq` - Raw sequencing reads (FASTQ/FASTA)

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--expected_telomere_length` | 8000 | Expected telomere length |
| `--phred_threshold` | 20 | Quality threshold for reads |
| `--length_threshold` | 0 | Minimum contig length |
| `--output_dir` | assemblmore_output | Output directory |
| `--skip_bam` | false | Skip BAM generation (PAF only) |

### Usage Examples

#### Basic assembly improvement:
```bash
./assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fastq
```

#### With custom parameters:
```bash
./assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fastq \
    --expected_telomere_length 10000 \
    --phred_threshold 25 \
    --output_dir my_results
```

#### Fast processing (skip BAM generation):
```bash
./assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fastq --skip_bam
```

#### Arguments in any order:
```bash
./assemblmore_pipeline.sh --output_dir results reference.fasta --phred_threshold 30 assembly.fasta reads.fastq
```

## Pipeline Steps

### 1. Reference Mapping
- Maps assembly contigs to reference genome using minimap2
- Determines optimal contig placement and orientation
- Filters contigs based on quality thresholds

### 2. Contig Placement
- Analyzes mapping results to determine contig order
- Handles complex regions and overlapping mappings
- Creates filtered contig placement file

### 3. Assembly Refinement  
- Generates improved assembly with corrected contig order
- Applies reverse complementation where needed
- Creates intermediate refined assembly

### 4. Gap Filling
- Maps reads to refined assembly
- Identifies and fills gaps with high-quality sequences
- Extends contigs using spanning reads

### 5. Quality Assessment
- Generates comprehensive assembly statistics
- Creates comparative analysis (before/after)
- Produces visualization plots and coverage analysis

## Output Files

The pipeline generates numerous output files in the specified output directory:

### Primary Outputs
- `assemblmore_final_assembly.fasta` - Final improved assembly
- `final_assembly_coverage.txt` - Coverage analysis
- `assembly_comparison.csv` - Comparative statistics

### Intermediate Files
- `*_mapped_contigs.paf` - Mapping alignments
- `filtered_*_contigs.tsv` - Filtered contig placements
- `ordered_and_oriented_*.fasta` - Intermediate assemblies
- `*.sam`/`*.bam` - Alignment files (if not skipped)

### Analysis and Visualization
- `*_nx_plot.png` - NX curve comparisons
- `*_length_distribution.png` - Contig length distributions
- `*_assembly_stats.csv` - Detailed statistics
- Various R-generated plots and analysis files

## Advanced Features

### Comparative Analysis

The pipeline automatically compares original and improved assemblies:

```bash
# Results include side-by-side statistics
./assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fastq
```

### Manual Comparative Analysis

```bash
# Compare multiple assemblies
Rscript assembly_stats.R \
    original.fasta:Original \
    improved.fasta:Improved \
    polished.fasta:Polished \
    output_directory
```

### Custom Assembly Statistics

```bash
# Generate detailed statistics for any assembly
Rscript assembly_stats.R assembly.fasta:MyAssembly output_dir
```

## Performance and Optimization

### Memory Usage
- Large assemblies may require substantial memory
- Consider using `--skip_bam` for faster processing
- Monitor system resources during execution

### Processing Time
- Runtime scales with assembly size and read count
- Typical runtime: 30 minutes to several hours
- Use multiple cores where possible (minimap2 threading)

### Quality Considerations
- Best results with high-quality reference genomes
- Optimal for assemblies with N50 ≥ 45kb
- Long-read data (ONT/PacBio) recommended

## Troubleshooting

### Common Issues

1. **Missing dependencies:**
   ```bash
   ./install_dependencies.sh
   ```

2. **Permission errors:**
   ```bash
   chmod +x *.sh
   ```

3. **Path issues:**
   ```bash
   # Use absolute paths for inputs
   ./assemblmore_pipeline.sh /full/path/to/reference.fasta ...
   ```

4. **Memory errors:**
   - Use smaller datasets for testing
   - Enable `--skip_bam` to reduce memory usage
   - Monitor system resources

### Debugging

Enable verbose mode for detailed logging:
```bash
./assemblmore_pipeline.sh -v reference.fasta assembly.fasta reads.fastq
```

Keep intermediate files for debugging:
```bash
./assemblmore_pipeline.sh -k reference.fasta assembly.fasta reads.fastq
```

## NOTE
All testing was done with relatively good datasets (N50 >= 45kb) on hermaphroditic or inbred worms. 

When adding /src/ to path, make sure to NOT use ~/ shell expansion. Use the $HOME value instead.

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/improvement`)
5. Create a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use AssemblMore in your research, please cite:

```
AssemblMore: A comprehensive genome assembly improvement pipeline
GitHub: https://github.com/TejusK123/assemblmore
```

## Support

- **Issues:** Report bugs and request features on [GitHub Issues](https://github.com/TejusK123/assemblmore/issues)
- **Documentation:** Check the `src/README.md` for detailed technical documentation
- **Examples:** See `example_comparative_analysis.md` for usage examples

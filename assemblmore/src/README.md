# AssemblMore - Genome Assembly Improvement Pipeline

AssemblMore is a comprehensive pipeline designed to improve genome assemblies created by next-generation assemblers like Flye and Canu. It uses reference-guided approaches and long-read data to enhance contig placement, orientation, and gap filling.

## Pipeline Scripts and Components

The AssemblMore pipeline consists of several interconnected scripts, each with a specific role in the assembly improvement process. The main orchestrator coordinates the execution of these individual components. Each component can be run individually.

### Core Pipeline Scripts

- **`assemblmore_pipeline.sh`** - Main orchestrator script that manages the entire workflow, coordinates execution of all pipeline steps with error handling, and provides verbose logging

- **`fill_gaps.sh`** - Wrapper for minimap2 sequence alignment that handles mapping between sequences, generates PAF alignment files and optionally BAM files, and manages output directories

- **`initial_assembly.sh`** - Creates initial refined assembly by coordinating contig placement analysis and reorienting contigs based on mapping results

- **`contig_placements.py`** - Analyzes PAF mappings to determine optimal contig placement, filters contigs based on mapping quality, and resolves contig orientations using network analysis

- **`span_contigs.sh`** - Simple bash wrapper that calls the Python gap filling implementation

- **`span_contigs.py`** - Core gap filling and assembly finalization script that performs sophisticated gap filling between contigs, implements telomere extension, and generates the final improved assembly

- **`assembly_stats.R`** - Comprehensive assembly quality analysis that calculates contiguity metrics, generates comparative plots, and implements quality assessment using clustering algorithms

### Utility and Supporting Scripts

- **`kmer_completeness.sh`** - K-mer based completeness analysis using Meryl and Merqury to provide quality metrics independent of reference genome

- **`install_dependencies.sh`** - Automated dependency installation for required bioinformatics tools, Python packages, and R libraries

- **`test_pipeline.sh`** - End-to-end pipeline testing script

- **`test_arguments.sh`** - Argument parsing validation

- **`test_argument_order.sh`** - Tests flexible argument ordering

### Pipeline Execution Flow

The pipeline follows this sequence: input validation → reference mapping → contig placement analysis → initial assembly creation → read mapping → gap filling → quality assessment and comparison

## Overview

The pipeline consists of four main steps:

1. **Reference Mapping**: Maps assembly contigs to a reference genome
2. **Contig Placement**: Determines optimal contig placement and orientation
3. **Assembly Refinement**: Creates an improved assembly with corrected contig order
4. **Gap Filling**: Maps reads to the refined assembly and fills remaining gaps

## Installation

### Prerequisites

Make sure you have the following tools installed:

- `minimap2` - For sequence alignment
- `samtools` - For SAM/BAM file manipulation
- `paftools.js` - For PAF format conversion (usually comes with minimap2)
- `python` - For running the analysis scripts
- `seqkit` - For FASTA/FASTQ manipulation

### Python Dependencies

The Python scripts require:
- `numpy`
- `pandas` 
- `biopython`
- `networkx`
- `more_itertools`
- `click`

Install with:
```bash
pip install numpy pandas biopython networkx more_itertools click
```

## Usage

### Basic Usage

```bash
./assemblmore_pipeline.sh reference_genome.fasta assembly.fasta nanopore_reads.fasta
```

### Advanced Usage

```bash
./assemblmore_pipeline.sh [OPTIONS] <reference_genome.fasta> <assembly.fasta> <nanopore_reads.fasta>
```

### Options

- `-o, --output-dir DIR`: Output directory (default: assemblmore_output)
- `-p, --preset1 PRESET`: Minimap2 preset for step 1 (default: asm20)
- `-P, --preset2 PRESET`: Minimap2 preset for step 3 (default: map-ont)
- `-N, --max-alignments N`: Maximum alignments per read (default: 5)
- `-t, --telomere-length N`: Expected telomere length for gap filling (default: 8000)
- `-l, --length-threshold N`: Minimum read length threshold (default: 0)
- `-q, --phred-threshold N`: Minimum Phred quality threshold (default: 20)
- `-k, --keep-intermediate`: Keep intermediate files for debugging
- `-v, --verbose`: Enable verbose output
- `-h, --help`: Show help message

### Examples

#### Basic assembly improvement:
```bash
./assemblmore_pipeline.sh reference.fasta my_assembly.fasta reads.fasta
```

#### With custom parameters for gap filling and quality control:
```bash
./assemblmore_pipeline.sh -t 10000 -q 30 -l 5000 reference.fasta my_assembly.fasta reads.fasta
```

#### Advanced example with all options:
```bash
./assemblmore_pipeline.sh \
    -o my_results \
    -p asm10 \
    -P map-hifi \
    -N 10 \
    -t 12000 \
    -l 1000 \
    -q 25 \
    -k -v \
    reference.fasta assembly.fasta reads.fasta
```

#### Keep intermediate files for debugging:
```bash
./assemblmore_pipeline.sh -k -v reference.fasta assembly.fasta reads.fasta
```

#### Arguments can be in any order:
```bash
# Options before positional arguments
./assemblmore_pipeline.sh -t 10000 -v reference.fasta assembly.fasta reads.fasta

# Options after positional arguments  
./assemblmore_pipeline.sh reference.fasta assembly.fasta reads.fasta -t 10000 -v

# Mixed order
./assemblmore_pipeline.sh -t 10000 reference.fasta -v assembly.fasta reads.fasta -q 30
```

### Running Individual Pipeline Components

For advanced users who need to run specific pipeline steps or customize the workflow:

#### Manual Step-by-Step Execution

```bash
# Step 1: Map assembly to reference
./fill_gaps.sh -ax asm20 reference.fasta assembly.fasta

# Step 2: Analyze contig placements and create initial assembly
./initial_assembly.sh assembly_to_reference.paf reads.fastq assembly.fasta

# Step 3: Map reads to refined assembly
./fill_gaps.sh -ax map-ont refined_assembly.fasta reads.fastq

# Step 4: Generate final assembly with gap filling
./span_contigs.sh read_alignments.paf contig_placements.tsv \
    refined_assembly.fasta reads.fastq \
    --expected_telomere_length 8000 --phred_threshold 20

# Step 5: Generate assembly statistics and comparisons
Rscript assembly_stats.R original_assembly.fasta final_assembly.fasta \
    --output-dir results/ --coverage-file coverage.txt
```

#### Individual Script Usage

**fill_gaps.sh** - Minimap2 wrapper with PAF output:
```bash
# Basic alignment
./fill_gaps.sh -ax asm20 reference.fasta query.fasta

# Skip BAM generation for speed
./fill_gaps.sh --no-bam -ax map-ont assembly.fasta reads.fastq

# Custom output directory
./fill_gaps.sh -o custom_output/ -ax asm20 ref.fasta assembly.fasta
```

**contig_placements.py** - Contig placement analysis:
```bash
# Basic usage
python contig_placements.py mappings.paf reads.fastq

# The script automatically generates placement TSV files
```

**span_contigs.py** - Gap filling and assembly finalization:
```bash
# Full parameter specification
python span_contigs.py alignments.paf placements.tsv assembly.fasta reads.fastq \
    --expected_telomere_length 10000 \
    --phred_threshold 25 \
    --length_threshold 1000 \
    --disable_telomere_extension  # For bacterial genomes

# Basic usage with defaults
python span_contigs.py alignments.paf placements.tsv assembly.fasta reads.fastq
```

**kmer_completeness.sh** - K-mer based quality assessment:
```bash
# Basic k-mer analysis
./kmer_completeness.sh reads.fastq assembly.fasta

# Custom k-mer size and output
./kmer_completeness.sh reads.fastq assembly.fasta --kmer-size 21 --output-dir kmer_results/
```

## Pipeline Steps in Detail

### Step 1: Reference Mapping
Maps the input assembly contigs to the reference genome using minimap2 with the specified preset (default: asm20).

**Input**: Reference genome, assembly contigs
**Output**: PAF alignment file

### Step 2: Contig Placement Analysis
Analyzes the mappings to determine:
- Which contigs belong to which chromosomes
- Optimal orientation for each contig
- Overlapping regions that need to be resolved

**Input**: PAF file, nanopore reads, assembly
**Output**: TSV file with contig placements, initial refined assembly

### Step 3: Read Mapping to Refined Assembly
Maps the original nanopore reads to the refined assembly to identify remaining gaps and misassemblies.

**Input**: Refined assembly, nanopore reads
**Output**: PAF alignment file

### Step 4: Final Assembly Generation
Performs final gap filling and assembly polishing based on the read mappings. This step uses several configurable parameters:

- **Expected Telomere Length** (`-t`): Expected length of telomeric regions for extension (default: 8000 bp)
- **Length Threshold** (`-l`): Minimum read length to consider for gap filling (default: 0 bp)
- **Phred Threshold** (`-q`): Minimum mapping quality score (Phred scale) for reads to be used (default: 20, equivalent to 99% accuracy)

**Input**: Read alignments, contig placements, refined assembly, reads
**Output**: Final improved assembly

## Parameter Guidelines

### Telomere Length (`-t`, `--telomere-length`)
- **Default**: 8000 bp
- **Purpose**: Sets the expected length of telomeric regions for gap filling and extension
- **Recommendation**: 
  - Use 8000-10000 for most eukaryotic genomes
  - Use 5000-8000 for smaller genomes or organisms with shorter telomeres
  - Use 10000-15000 for genomes known to have very long telomeres

### Length Threshold (`-l`, `--length-threshold`)
- **Default**: 0 bp (no filtering)
- **Purpose**: Filters out reads shorter than this threshold for gap filling
- **Recommendation**:
  - Use 0 to include all reads (default)
  - Use 1000-5000 to focus on longer, more reliable reads
  - Higher values improve accuracy but may reduce coverage in gap regions

### Phred Threshold (`-q`, `--phred-threshold`)
- **Default**: 20 (99% accuracy)
- **Purpose**: Minimum mapping quality score for reads to be considered reliable
- **Recommendation**:
  - Use 10-15 for noisy data or when coverage is limited
  - Use 20-30 for standard quality data (recommended)
  - Use 30+ for high-quality data where precision is critical
  - Higher values improve accuracy but may exclude useful reads

## Output Files

The pipeline generates several output files in the specified output directory:

- `assemblmore_final_assembly.fasta` - The final improved assembly
- `filtered_by_<reference>_contigs.tsv` - Contig placement information
- `ordered_and_oriented_to_<reference>_assembly.fasta` - Initial refined assembly
- Various intermediate alignment files (if `-k` option is used)

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Make sure all required tools are installed and in your PATH
2. **File format issues**: Ensure input files are in the correct format (FASTA for genomes/assemblies)
3. **Memory issues**: Large assemblies may require significant RAM for processing
4. **Path issues**: Use absolute paths if you encounter file not found errors

### Debugging

Use the `-v` (verbose) and `-k` (keep intermediate files) options to get more detailed output and preserve intermediate files for inspection.

### Log Files

The pipeline provides timestamped logging to help track progress and identify issues.

## Performance Notes

- Processing time depends on assembly size, read count, and system resources
- Large assemblies (>1GB) may take several hours to process
- Consider using a system with adequate RAM (16GB+ recommended for large genomes)

## Contact

For issues or questions about AssemblMore, please check the repository issues or create a new issue with detailed information about your problem.

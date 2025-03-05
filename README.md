# Genomic TE Simulation Toolkit

A comprehensive toolkit for simulating the insertion, amplification, and evolution of transposable elements (TEs) in genomic sequences.

## Overview

This toolkit provides a collection of Python scripts for computational genomics research focused on transposable elements. It enables researchers to simulate how TEs insert into genomes, amplify over time, and evolve through various mutation events, producing realistic genomic sequences with annotated TE content for downstream analysis.

## Features

- Generate structured random genomes with realistic GC content distribution and chromosome features
- Simulate TE insertions with preferences for specific genomic regions or GC content
- Model genome evolution through SNPs and structural variations (deletions, insertions, inversions)
- Analyze TE sequence similarity and estimate evolutionary time
- Generate comprehensive statistics and visualizations of TE distribution
- Parallel processing support for handling large genomes efficiently
- Multi-stage evolutionary simulation with configurable parameters

## Requirements

- Python 3.6 or higher
- Required Python packages:
  - NumPy
  - Matplotlib
  - BioPython
  - Pandas (for some visualization scripts)
- Optional: R (for advanced visualizations)

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/username/genomic-te-simulation.git
   cd genomic-te-simulation
   ```

2. Install required Python packages:
   ```
   conda create -f env/environment.yml
   ```

3. Make the main script executable:
   ```
   chmod +x run_simulation.sh
   ```

## Usage

### Quick Start

Run a complete simulation with default parameters:

```bash
./run_simulation.sh -t path/to/te_library.fa -o output_directory
```

### Main Workflow

The complete workflow consists of four main steps:

1. **Generate a base genome** with realistic chromosome structure
2. **Insert TEs** into the genome according to specified preferences
3. **Simulate evolution** by introducing mutations and structural variations
4. **Analyze the results** to generate annotations, statistics, and visualizations

### Command-line Options

```bash
./run_simulation.sh [options]
```

Main options:
- `-l, --length GENOME_LENGTH`: Initial genome length (default: 50,000,000)
- `-p, --percentage TE_PERCENTAGE`: Target TE percentage (default: 30)
- `-e, --evolution EVOLUTION_LEVEL`: Evolution level (low/medium/high, default: medium)
- `-s, --stages STAGES`: Number of evolutionary stages (default: 3)
- `-t, --te_library TE_LIBRARY`: TE library file (required)
- `-o, --output OUTPUT_DIR`: Output directory (default: te_simulation_output)
- `-d, --distribution TE_DISTRIBUTION`: TE insertion stage distribution (comma-separated percentages, default: 60,30,10)
- `--seed SEED`: Random seed (default: current timestamp)
- `--gc GC_CONTENT`: Base GC content (default: 0.5)
- `--gc_gradient`: Enable GC content gradient
- `--mutation_rate MUTATION_RATE`: Per base per generation mutation rate (default: 1e-8)
- `--num_processes NUM_PROCESSES`: Number of processes for parallel processing (default: automatic)
- `--fast`: Use fast TE insertion algorithm (ignores region preferences and GC content)

For full options list:
```bash
./run_simulation.sh --help
```

### Using Individual Scripts

You can also run individual components of the workflow:

1. **Generate a base genome**:
   ```bash
   python generate_base_genome.py -l 50000000 -o base_genome.fa --gc 0.5
   ```

2. **Insert TEs**:
   ```bash
   python optimized_te_insertion.py -g base_genome.fa -t te_library.fa -p 30 -o te_inserted
   ```

3. **Simulate evolution**:
   ```bash
   python simulate_evolution.py -g te_inserted_genome.fa -a te_inserted_te_annotations.bed -l medium -o evolved
   ```

4. **Generate annotations and statistics**:
   ```bash
   python generate_annotation.py -g evolved_genome.fa -a evolved_annotations.bed -o te_annotation
   ```

5. **Analyze TE similarity**:
   ```bash
   python te_similarity.py -t te_library.fa -s evolved_te_sequences.fa -a evolved_annotations.bed -o te_similarity
   ```

## Script Descriptions

- `generate_base_genome.py`: Generates an initial random genome sequence with specific GC content and chromosome structures
- `optimized_te_insertion.py`: Inserts TEs into a genome with optimized performance for large genomes
- `simulate_evolution.py`: Introduces SNPs and structural variations to simulate genomic evolution
- `generate_annotation.py`: Generates annotation files and statistical reports
- `te_similarity.py`: Analyzes TE sequence similarity and estimates evolutionary time
- `time_evolution.py`: Simulates multiple evolutionary stages including TE insertion and mutation
- `genomic_partition.py`: Provides parallel genomic partition processing
- `parallel_processing.py`: Provides reusable parallel processing functions
- `parameter_validation.py`: Validates input parameters for different scripts
- `visualization_tools.py`: Generates visualizations of chromosome structure and TE distribution
- `run_simulation.sh`: Main script to run the complete simulation workflow

## Output Files

For each simulation, the toolkit produces:

- **Genome sequences** (FASTA format)
- **TE annotations** (BED and GFF3 formats)
- **TE sequences** extracted from the genome
- **Statistical reports** on TE content, distribution, and evolution
- **Visualizations** of TE distribution, similarity, and evolutionary time
- **Evolution logs** showing mutation events

## Example

Here's a basic example of running a simulation with a custom TE library:

```bash
./run_simulation.sh \
  -t example_data/te_library.fa \
  -l 10000000 \
  -p 25 \
  -e medium \
  -s 2 \
  -o example_output \
  --gc 0.45 \
  --gc_gradient \
  --seed 12345
```

This will:
1. Generate a 10Mb random genome with 45% GC content and a gradient
2. Perform 2 stages of TE insertion, targeting 25% overall TE content
3. Apply medium-level evolutionary changes
4. Save all results to the `example_output` directory

## Citation

If you use this toolkit in your research, please cite:

PanTE: A Comprehensive Framework for Transposable Element Discovery in Graph-based Pangenomes[(https://www.researchgate.net/profile/Cao-Shuo-7/publication/388740613_PanTE_A_Comprehensive_Framework_for_Transposable_Element_Discovery_in_Graph-based_Pangenomes/links/67aa0438207c0c20fa836d65/PanTE-A-Comprehensive-Framework-for-Transposable-Element-Discovery-in-Graph-based-Pangenomes.pdf)]

## License

This project is licensed under the [MIT License](LICENSE)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

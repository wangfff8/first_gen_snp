# VSAT: Virus SNP Analysis Tool

**VSAT (Virus SNP Analysis Tool)** is a simple and efficient Python package for studying Single Nucleotide Polymorphisms (SNPs) in viral genomes using first-generation sequencing data.

## Features

- Automates the initial assembly of raw sequencing data using CAP3.
- Identifies SNPs by comparing the assembled genome to a reference sequence.
- Provides detailed information about mutations, including their impact on amino acid coding.
- Designed as an importable library, allowing for easy integration into custom workflows.

## Requirements

- Python 3.10+
- [uv](https://github.com/astral-sh/uv) (recommended for environment management)
- **[CAP3](https://faculty.sites.iastate.edu/xqhuang/cap3-assembly-program)**: Must be installed and available in your system's `PATH`.
- **[SeqKit](https://github.com/shenwei356/seqkit/releases)**: Must be installed and available in your system's `PATH`.

## Installation

1. **Clone the repository:**

    ```bash
    git clone <repository_url>
    cd first_gen_snp
    ```

2. **Create and activate a virtual environment:**

    ```bash
    uv venv
    # On Windows (Git Bash)
    source .venv/Scripts/activate
    # On macOS/Linux
    source .venv/bin/activate
    ```

3. **Install the package in editable mode:**
    This makes the `vsat` package available for import in your Python environment.

    ```bash
    uv pip install -e .
    ```

## Project Structure

The project uses a modern `src` layout:

- `src/vsat/`: Main source code for the VSAT package.
- `data/`: Contains reference genome data.
  - `data/genomes/`: Stores reference genomes (.fasta) and annotations (.gff3).
    - `data/ref_genome.json`: Maps locus names to available genome files.
- `scripts/`: Contains example shell scripts for running the workflows.
- `pyproject.toml`: Defines project metadata and dependencies.

## Usage

The VSAT package can be used either directly from the command line or as a library in your own Python scripts.

### Command-Line Usage

#### Step 1: Data Preparation

1. **Reference Genomes**:
    - Place your reference virus genomes (`.fasta` files) and their corresponding annotation files (`.gff3` files) into the `data/genomes/` directory.
    - Update the `data/ref_genome.json` file to map a short name (locus) to your new genome files.

2. **Raw Data**:
    - Organize your raw sequencing files (`.seq`, `.ab1`, etc.) in a dedicated directory (e.g., `/path/to/raw_data`).

3. **ID Map File**:
    - Create an ID map file (e.g., `id_map.xls`). This is a two-column, tab-separated file that maps a unique identifier from your raw data filenames to a desired sample name.

    ```plaintext
    # Column 1: Identifier in filename   Column 2: Desired sample name
    320240824010  Me2024074S05_01-3
    ```

#### Step 2: Sequence Assembly

The genome assembly workflow is run via the `vsat.genome_assembler` module.

- A template script is provided in `scripts/run_genome_assembler.sh`. You can copy, modify, and run this script.
- **Important**: Before running, make sure your virtual environment is activated.

Example command:

```bash
python -m vsat.genome_assembler \
    --id_map /path/to/id_map.xls \
    --raw_data /path/to/raw_data \
    --split_data /path/to/split_data \
    --assembly_dir /path/to/assembly_dir \
    
```

#### Step 3: Manual Finishing

After the initial assembly, you will find the assembled sequences in the output directory (`assemble/*/*.cap.contigs.split`).

These contigs require manual alignment and finishing using tools like DNAMAN, especially by comparing them against the raw `.ab1` trace files and the reference genome to produce a complete, high-quality viral genome.

#### Step 4: SNP Detection

Once you have a finished assembly (saved as a `.seq` file in a directory), use the `vsat.variant_caller` module to detect SNPs against the reference genome.

To see a list of available reference genomes (defined in `data/ref_genome.json`):

```bash
python -m vsat.variant_caller --list
```

To run SNP detection:

```bash
# -g: The locus name from the reference genome list
# -s: The path to the directory containing your final assembled .seq files
python -m vsat.variant_caller -g <locus_name> -s /path/to/final_assembly_dir
```

The results, including `snp_result.txt`, `snp_result.xls`, and alignment statistics, will be generated in your current working directory.

### Library Usage

You can also import and use the workflow functions directly in your Python scripts.

#### Example: Genome Assembly

```python
from vsat import genome_assembler

# Define paths for the assembly workflow
id_map = "/path/to/id_map.xls"
raw_data = "/path/to/raw_data"
split_dir = "/path/to/split_data"
assembly_dir = "/path/to/assembly_dir"

# Run the assembly process
genome_assembler.run_genome_assembly(
    id_map_file=id_map,
    raw_data_dir=raw_data,
    split_data_dir=split_dir,
    assembly_dir=assembly_dir,
)
```

#### Example: SNP Analysis

```python
from vsat import variant_caller

# Define paths for the SNP analysis
assembled_dir = "/path/to/final_assembly_dir"
ref_fasta = "data/genomes/<LOCUS_NAME>.fasta"
ref_gff = "data/genomes/<LOCUS_NAME>.gff3"

# Run the SNP analysis
variant_caller.run_snp_analysis(
    ref_genome_fasta_file=ref_fasta,
    ref_genome_gff_file=ref_gff,
    assembled_sequences_dir=assembled_dir,
)
```

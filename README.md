# VSAT: Virus SNP Analysis Tool

**VSAT (Virus SNP Analysis Tool)** is a simple and efficient Python package for
studying Single Nucleotide Polymorphisms (SNPs) in viral genomes using
first-generation sequencing data. It provides an end-to-end workflow from raw
data to a final, interactive HTML report.

## Features

- Automates the initial assembly of raw sequencing data using **CAP3**.
- Identifies SNPs by comparing the assembled genome to a reference sequence.
- Performs multiple sequence alignment for genomes, genes, and proteins using
  **MAFFT**.
- Generates a comprehensive, interactive **HTML report** visualizing SNP results
  and sequence alignments.
- Provides detailed information about mutations, including their impact on amino
  acid coding.
- Designed as an importable library, allowing for easy integration into custom
  workflows.

## Requirements

- Python 3.12+
- [uv](https://github.com/astral-sh/uv) (recommended for environment management)
- **[CAP3](https://faculty.sites.iastate.edu/xqhuang/cap3-assembly-program)**:
  Must be installed and available in your system's `PATH`.
- **[SeqKit](https://github.com/shenwei356/seqkit/releases)**: Must be installed
  and available in your system's `PATH`.
- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)**: Must be installed and
  available in your system's `PATH`.

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

3. **Install the package in editable mode:** This makes the `vsat` package
    available for import in your Python environment.

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

The VSAT package can be used either directly from the command line or as a
library in your own Python scripts.

### Command-Line Usage

#### Step 1: Data Preparation

1. **Reference Genomes**:
    - Place your reference virus genomes (`.fasta` files) and their
      corresponding annotation files (`.gff3` files) into the `data/genomes/`
      directory.
    - Update the `data/ref_genome.json` file to map a short name (locus) to your
      new genome files.

2. **Raw Data**:
    - Organize your raw sequencing files (`.seq`, `.ab1`, etc.) in a dedicated
      directory (e.g., `/path/to/raw_data`).

3. **ID Map File**:
    - Create an ID map file (e.g., `id_map.xls`). This is a two-column,
      tab-separated file that maps a unique identifier from your raw data
      filenames to a desired sample name.

    ```plaintext
    # Column 1: Identifier in filename   Column 2: Desired sample name
    32024082401051  Me2024074S05_01-3
    ```

#### Step 2: Sequence Assembly

The genome assembly workflow is run via the `vsat.genome_assembler` module. A
template script is provided in `scripts/genome_assembly.sh`.

Example command:

```bash
python -m vsat.genome_assembler \
    --id_map /path/to/id_map.xls \
    --raw_data /path/to/raw_data \
    --split_data /path/to/split_data \
    --assembly_dir /path/to/assembly_dir
```

After this step, you may need to manually inspect and finish the assembled
contigs (found in `assembly_dir`) using tools like DNAMAN to produce a
high-quality viral genome sequence (`.seq` file).

#### Step 3: SNP Detection and Reporting

Once you have a finished assembly (saved as a `.seq` file in a directory), use
the `vsat.variant_caller` module to detect SNPs and generate a full analysis
report. A template script is provided in `scripts/snp_detection.sh`.

First, list available reference genomes:

```bash
python -m vsat.variant_caller --list
```

Then, run the SNP detection workflow:

```bash
# -g: The locus name from the reference genome list
# -s: The path to the directory containing your final assembled .seq files
# -o: The directory where all results will be stored
python -m vsat.variant_caller \
    -g <locus_name> \
    -s /path/to/final_assembly_dir \
    -o /path/to/output_results_dir
```

The output directory will contain:

- `viral_snp_report.html`: The main interactive report.
- `snp_result.txt` & `snp_result.xls`: Detailed SNP data.
- `alignment_statistics.xls`: Sequence identity metrics.
- `processed_sequences/`: A structured directory with all intermediate sequences
  and MAFFT alignment results for genomes, genes, and proteins.

### Library Usage

You can also import and use the workflow functions directly in your Python
scripts.

#### Example: Genome Assembly

```python
from vsat import genome_assembler

genome_assembler.run_genome_assembly(
    id_map_file="/path/to/id_map.xls",
    raw_data_dir="/path/to/raw_data",
    split_data_dir="/path/to/split_data",
    assembly_dir="/path/to/assembly_dir",
)
```

#### Example: SNP Analysis

```python
from pathlib import Path
from vsat import variant_caller

# Define paths for the SNP analysis
ref_fasta = Path("data/genomes/<LOCUS_NAME>.fasta")
ref_gff = Path("data/genomes/<LOCUS_NAME>.gff3")
assembled_dir = Path("/path/to/final_assembly_dir")
output_dir = Path("/path/to/output_results_dir")

# Run the SNP analysis
variant_caller.run_snp_analysis(
    ref_genome_fasta_file=ref_fasta,
    ref_genome_gff_file=ref_gff,
    assembled_sequences_dir=assembled_dir,
    output_dir=output_dir,
)
```

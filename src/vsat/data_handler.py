"""
This module provides utility functions for data handling and file I/O.

It includes functions for creating directories, parsing various bioinformatics
file formats (FASTA, GFF3, custom ID maps), and organizing raw sequencing
data into a structured format for analysis.
"""

import re
import json
import shutil
from itertools import chain
from pathlib import Path


def mkdir(dir_path: str | Path) -> None:
    """Create a directory, including any necessary parent directories."""
    Path(dir_path).mkdir(parents=True, exist_ok=True)


def load_reference_genomes(genome_file: str | Path) -> dict[str, list[str]]:
    """
    Load reference genome metadata from a JSON file.

    Args:
        genome_file: The path to the JSON file containing genome metadata.

    Returns:
        A dictionary parsed from the JSON file.
    """
    with open(genome_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def read_genome_sequence(seq_file: str | Path) -> tuple[str, str]:
    """
    Read a single FASTA-like file and standardize its format in place.

    This function reads a sequence file, extracts the first word of the header
    as the sequence ID, and consolidates the sequence onto a single line.

    Note:
        This function has a significant side effect: it overwrites the original
        file (`seq_file`) with the standardized FASTA content.

    Args:
        seq_file: The path to the sequence file.

    Returns:
        A tuple containing the sequence ID and the full sequence in uppercase.
    """
    seq_file_path = Path(seq_file)

    with open(seq_file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    seq_id = ""
    seq_parts = []
    for line in lines:
        if line.startswith(">"):
            seq_id = line.strip().split()[0][1:]
        else:
            seq_parts.append(line.strip().upper())

    # If no FASTA header is found, use the file's stem as the ID.
    if not seq_id:
        seq_id = seq_file_path.stem

    sequence = "".join(seq_parts)

    # Overwrite the original file with the standardized format.
    with open(seq_file_path, "w", encoding="utf-8") as f:
        f.write(f">{seq_id}\n{sequence}\n")

    return seq_id, sequence


def read_assembled_sequences(directory: str | Path) -> dict[str, str]:
    """
    Read all .seq files in a directory and load them into a dictionary.

    Args:
        directory: The path to the directory containing the .seq files.

    Returns:
        A dictionary mapping sequence names to their corresponding sequences.
    """
    genome_seq_dict = {}
    # Sort the glob results to ensure a consistent processing order.
    for genome_file in sorted(Path(directory).glob("*.seq")):
        seq_name, seq = read_genome_sequence(seq_file=genome_file)
        genome_seq_dict[seq_name] = seq
    return genome_seq_dict


def parse_gff(gff3_file: str | Path) -> dict[str, tuple[int, int]]:
    """
    Parse a GFF3 file to extract coding regions (CDS) and their product names.

    Args:
        gff3_file: The path to the GFF3 annotation file.

    Returns:
        A dictionary mapping product names to a tuple of (start, end) coordinates.
    """
    coding_dict = {}
    with open(gff3_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue

            # We are only interested in 'CDS' features.
            if fields[2] == "CDS":
                # Extract the product name from the attributes field.
                if match := re.search('product=(.*?)(?:;|$)', fields[8]):
                    product_name = match.group(1)
                    if product_name not in coding_dict:
                        coding_dict[product_name] = (int(fields[3]), int(fields[4]))

    return coding_dict


def read_id_map(id_map_file: str | Path) -> dict[str, str]:
    """
    Read a two-column, space-separated ID map file.

    Args:
        id_map_file: The path to the ID map file.

    Returns:
        A dictionary mapping source IDs from the first column to target sample
        names from the second column.
    """
    id_map_dict = {}
    if Path(id_map_file).is_file():
        with open(id_map_file, 'r', encoding='utf-8') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 2:
                    id_map_dict[fields[0]] = fields[1]
    return id_map_dict


def split_data(id_map: dict[str, str], raw_data_path: str | Path, out_path: str | Path) -> None:
    """
    Organize raw sequencing data into sample-specific directories based on an ID map.

    This function walks through the `raw_data_path`, finds all relevant files
    (.seq, .ab1, .pdf), and copies them into subdirectories under `out_path`,
    named according to the sample names in the `id_map`.

    Note:
        The script will print an error and exit if the number of .seq and .ab1
        files in the source directory do not match.

    Args:
        id_map: A dictionary mapping file identifiers to sample names.
        raw_data_path: The path to the directory containing raw sequencing files.
        out_path: The base directory for the organized, sample-specific output.
    """
    raw_data_path = Path(raw_data_path)
    out_path = Path(out_path)
    all_files = list(raw_data_path.rglob("*"))
    seq_list = [f for f in all_files if f.suffix == '.seq']
    ab1_list = [f for f in all_files if f.suffix == '.ab1']
    pdf_list = [f for f in all_files if f.suffix == '.pdf']

    if not seq_list or not ab1_list or len(seq_list) != len(ab1_list):
        print(
            'Error: Input directory must contain matching numbers of .seq and '
            '.ab1 files.'
        )
        exit(1)

    for file_id, sample_name in id_map.items():
        sample_out_path = out_path / sample_name
        # Create subdirectories for each file type within the sample folder.
        seq_out_path = sample_out_path / 'seq'
        ab1_out_path = sample_out_path / 'ab1'
        mkdir(dir_path=seq_out_path)
        mkdir(dir_path=ab1_out_path)

        pdf_out_path = None
        if pdf_list:
            pdf_out_path = sample_out_path / 'pdf'
            mkdir(dir_path=pdf_out_path)

        # Copy files that match the current file_id into the corresponding folders.
        for file_ in chain(seq_list, ab1_list, pdf_list):
            if re.search(file_id, file_.name):
                if file_.suffix == '.seq':
                    shutil.copy(file_, seq_out_path)
                elif file_.suffix == '.ab1':
                    shutil.copy(file_, ab1_out_path)
                elif file_.suffix == '.pdf' and pdf_out_path:
                    shutil.copy(file_, pdf_out_path)


def combine_seq2fasta(split_data_path: str | Path, assembly_dir: str | Path) -> None:
    """
    Combine all .seq files for a sample into a single multi-FASTA file.

    This function serves as a preprocessing step for assemblers like CAP3, which
    expect a single input file containing all reads for a given sample.

    Args:
        split_data_path: The directory containing sample-specific subdirectories,
                         each with its own `seq` folder of .seq files.
        assembly_dir: The base directory where the output multi-FASTA files will
                      be stored, inside sample-specific subfolders.
    """
    split_data_path = Path(split_data_path)
    assembly_dir = Path(assembly_dir)
    for sample_path in split_data_path.iterdir():
        if not sample_path.is_dir():
            continue
        sample_name = sample_path.name
        sample_seq_dir = sample_path / 'seq'
        if not sample_seq_dir.is_dir():
            continue

        # Create the output directory for the sample in the assembly folder.
        sample_out_path = assembly_dir / sample_name
        mkdir(dir_path=sample_out_path)
        output_fasta_file = sample_out_path / f"{sample_name}.fasta"
        with open(output_fasta_file, "w", encoding="utf-8") as o:
            seq_files = sorted(sample_seq_dir.glob("*.seq"))
            for seq_file in seq_files:
                seq_name = seq_file.stem
                with open(seq_file, "r", encoding="utf-8") as r:
                    o.write(f">{seq_name}\n")
                    o.write(r.read().strip() + "\n")


def convert_seq2fasta(in_path: str | Path, out_path: str | Path) -> None:
    """
    Convert individual .seq files into properly formatted .fasta files.

    This function traverses a directory of split sample data, finds all .seq
    files, and creates a corresponding .fasta file for each one in the output
    directory, adding a FASTA header based on the original filename.

    Args:
        in_path: The root directory containing the split sample data (e.g., 'split_data').
        out_path: The root directory where new sample directories and .fasta
                  files will be created.
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    for seq_file_fullname in in_path.rglob("*.seq"):
        # Assumes a directory structure of .../sample_name/seq/
        sample_name = seq_file_fullname.parent.parent.name
        sample_out_path = out_path / sample_name
        mkdir(dir_path=sample_out_path)

        seq_name = seq_file_fullname.stem
        out_seq_file = sample_out_path / f"{seq_name}.fasta"

        with (
            open(seq_file_fullname, "r", encoding="utf-8") as r,
            open(out_seq_file, "w", encoding="utf-8") as o,
        ):
            o.write(f">{seq_name}\n")
            o.write(r.read())

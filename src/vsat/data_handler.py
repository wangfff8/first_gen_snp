#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import re
import json
import shutil
from itertools import chain
from pathlib import Path


def mkdir(dir_path: str | Path) -> None:
    """Creates a directory if it does not already exist."""
    Path(dir_path).mkdir(parents=True, exist_ok=True)


def load_reference_genomes(genome_file: str | Path) -> dict[str, list[str]]:
    """Loads reference genome metadata from a JSON file."""
    with open(genome_file, 'r', encoding='utf-8') as f:
        return json.load(f)


def read_genome_sequence(seq_file: str | Path) -> tuple[str, str]:
    """
    Reads a single sequence file in FASTA-like format.

    Args:
        seq_file: The path to the sequence file.

    Returns:
        A tuple containing the sequence ID and the sequence content in uppercase.
    """
    with open(seq_file, "r", encoding="utf-8") as f:
        seq_id = ""
        seq = []
        for line in f:
            if line.startswith(">"):
                seq_id = line.strip().split()[0][1:]
            else:
                seq.append(line.strip().upper())
        
    return seq_id, ''.join(seq)


def read_assembled_sequences(directory: str | Path) -> dict[str, str]:
    """
    Reads all .seq files in a given directory into a dictionary.

    Args:
        directory: The path to the directory containing .seq files.

    Returns:
        A dictionary mapping sequence names to their corresponding sequences.
    """
    genome_seq_dict = {}
    for genome_file in Path(directory).glob("*.seq"):
        seq_name, seq = read_genome_sequence(genome_file)
        genome_seq_dict[seq_name] = seq
    return genome_seq_dict


def parse_gff(gff3_file: str | Path) -> dict[str, tuple[int, int]]:
    """
    Parses a GFF3 file and extracts product names and their coding regions (CDS).

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

            if fields[2] == "CDS":
                if match := re.search('product=(.*?)(?:;|$)', fields[8]):
                    product_name = match.group(1)
                    if product_name not in coding_dict:
                        coding_dict[product_name] = (int(fields[3]), int(fields[4]))

    return coding_dict


def read_id_map(id_map_file: str | Path) -> dict[str, str]:
    """
    Reads a two-column (source_id, target_name) ID map file.

    Args:
        id_map_file: The path to the ID map file.

    Returns:
        A dictionary mapping source IDs to target sample names.
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
    Splits and organizes raw sequencing data based on an ID map.

    Walks through the raw_data_path, finds all .seq, .ab1, and .pdf files,
    and copies them into sample-specific subdirectories in the out_path.

    Args:
        id_map: A dictionary mapping file identifiers to sample names.
        raw_data_path: The path to the directory with raw sequencing files.
        out_path: The base directory for the organized output.
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
        # Create subdirectories for each file type
        seq_out_path = sample_out_path / 'seq'
        ab1_out_path = sample_out_path / 'ab1'
        mkdir(seq_out_path)
        mkdir(ab1_out_path)

        pdf_out_path = None
        if pdf_list:
            pdf_out_path = sample_out_path / 'pdf'
            mkdir(pdf_out_path)

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
    Combines individual .seq files for each sample into a single FASTA file.

    This is a necessary preprocessing step for the CAP3 assembler.

    Args:
        split_data_path: The directory where split_data organized the files.
        assembly_dir: The base directory where the final FASTA files will be written.
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

        # Create the output directory for the sample in the assembly folder
        sample_out_path = assembly_dir / sample_name
        mkdir(sample_out_path)
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
    Converts individual .seq files for each sample into separate FASTA files.

    This function traverses a directory structure organized by sample, finds all
    .seq files, and creates a corresponding .fasta file for each one in the
    output directory.

    Args:
        in_path: The root directory containing the split sample data.
        out_path: The root directory where new sample directories and .fasta
            files will be created.
    """
    in_path = Path(in_path)
    out_path = Path(out_path)
    for seq_file_fullname in in_path.rglob("*.seq"):
        # Assumes a directory structure of .../sample_name/seq/
        sample_name = seq_file_fullname.parent.parent.name
        sample_out_path = out_path / sample_name
        mkdir(sample_out_path)

        seq_name = seq_file_fullname.stem
        out_seq_file = sample_out_path / f"{seq_name}.fasta"

        with (
            open(seq_file_fullname, "r", encoding="utf-8") as r,
            open(out_seq_file, "w", encoding="utf-8") as o,
        ):
            o.write(f">{seq_name}\n")
            o.write(r.read())

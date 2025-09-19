#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import argparse
import subprocess
from pathlib import Path
from . import data_handler


def _run_cap3_for_samples(assemble_dir: str | Path, sample_list: list[str]) -> None:
    """
    Runs the core assembly commands (CAP3 and SeqKit) for each sample.

    This function iterates through a list of samples, runs CAP3 to assemble
    contigs, and then uses SeqKit to split the resulting contigs into
    separate files. It requires CAP3 and SeqKit to be in the system's PATH.

    Args:
        assemble_dir: The base directory where sample subdirectories are located.
        sample_list: A list of sample names to process.
    """
    print("Starting assembly process...")
    assemble_dir = Path(assemble_dir)
    for sample in sample_list:
        sample_dir = assemble_dir / sample
        data_handler.mkdir(sample_dir)

        fasta_file = sample_dir / f"{sample}.fasta"
        contigs_file = fasta_file.with_suffix(".fasta.cap.contigs")

        print(f"  Running CAP3 for sample: {sample}... and splitting...")
        subprocess.run(["cap3", fasta_file.name], cwd=sample_dir, check=True)
        subprocess.run(
            ["seqkit", "split", "-i", contigs_file.name], cwd=sample_dir, check=True
        )

    print("Assembly process finished.")


def run_genome_assembly(
    raw_data_dir: str | Path,
    id_map_file: str | Path,
    split_data_dir: str | Path,
    assemble_dir: str | Path,
) -> None:
    """
    Orchestrates the entire genome assembly workflow.

    This workflow consists of four main steps:
    1. Loading the sample ID map.
    2. Splitting the raw sequencing data into sample-specific directories.
    3. Combining individual sequence files into a single FASTA file per sample.
    4. Running the CAP3 assembler and splitting the resulting contigs.

    Args:
        raw_data_dir: Directory containing the raw sequencing files (.seq, .ab1).
        id_map_file: Path to the ID map file mapping file IDs to sample names.
        split_data_dir: Directory to store the intermediate split data.
        assemble_dir: Directory to store the final assembled contigs.
    """
    print("Loading ID map...")
    id_map_dict = data_handler.read_id_map(id_map_file)
    sample_list = list(id_map_dict.values())
    print(f"Found {len(sample_list)} samples.")

    print("Splitting raw data...")
    data_handler.split_data(id_map_dict, raw_data_dir, split_data_dir)

    print("Combining sequences into FASTA files...")
    data_handler.combine_seq2fasta(split_data_dir, assemble_dir)

    _run_cap3_for_samples(assemble_dir, sample_list)


def main():
    """
    Parses command-line arguments and executes the genome assembly workflow.

    This function serves as the command-line entry point for the script.
    It parses arguments for input/output directories and the ID map,
    then calls the main assembly workflow function.
    """
    parser = argparse.ArgumentParser(
        description="Assemble raw sequencing data using CAP3."
    )
    parser.add_argument(
        "--rawdata", required=True, help="Directory containing raw sequencing files."
    )
    parser.add_argument(
        "--split_data",
        required=True,
        help="Directory to store intermediate split data.",
    )
    parser.add_argument(
        "--assemble_dir", required=True, help="Directory to store assembled contigs."
    )
    parser.add_argument("--id_map", required=True, help="Path to the ID map file.")
    args = parser.parse_args()

    run_genome_assembly(
        raw_data_dir=args.rawdata,
        id_map_file=args.id_map,
        split_data_dir=args.split_data,
        assemble_dir=args.assemble_dir,
    )


if __name__ == "__main__":
    main()

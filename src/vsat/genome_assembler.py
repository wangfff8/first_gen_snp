#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
This module orchestrates the genome assembly workflow for raw sequencing data.

It automates the process of organizing raw data, running the CAP3 assembler,
and preparing the resulting contigs for further analysis. This script can be
executed directly from the command line.

The workflow relies on two external command-line tools that must be installed
and available in the system's PATH:
- cap3: For sequence assembly.
- seqkit: For splitting contig files.
"""

import argparse
import subprocess
from pathlib import Path
from . import data_handler, logger as vsat_logger

# Set up logging with colored output.
logger = vsat_logger.setup_logger(__name__)


def _run_cap3_for_samples(assembly_dir: str | Path, sample_list: list[str]) -> None:
    """
    Run the CAP3 assembler and SeqKit splitter for each sample.

    This function iterates through a list of sample names, creates a dedicated
    directory for each, and then executes the core assembly commands: `cap3`
    to assemble contigs and `seqkit split` to separate them into individual
    files.

    Args:
        assembly_dir: The base directory where sample-specific assembly folders
                      will be created and results will be stored.
        sample_list: A list of sample names to be processed.
    """
    logger.info("Starting assembly process...")
    assembly_dir = Path(assembly_dir)
    for sample in sample_list:
        sample_dir = assembly_dir / sample
        data_handler.mkdir(sample_dir)

        fasta_file = sample_dir / f"{sample}.fasta"
        contigs_file = fasta_file.with_suffix(".fasta.cap.contigs")

        logger.info(f"  Running CAP3 for sample: {sample}... and splitting...")
        subprocess.run(["cap3", fasta_file.name], cwd=sample_dir, check=True)
        subprocess.run(
            ["seqkit", "split", "-i", contigs_file.name], cwd=sample_dir, check=True
        )

    logger.info("Assembly process finished.")


def run_genome_assembly(
    id_map_file: str | Path,
    raw_data_dir: str | Path,
    split_data_dir: str | Path,
    assembly_dir: str | Path,
) -> None:
    """
    Orchestrate the end-to-end genome assembly workflow.

    This function manages the entire assembly pipeline, which includes:
    1. Reading the sample ID map to associate file IDs with sample names.
    2. Splitting raw sequencing data into sample-specific subdirectories.
    3. Consolidating individual sequence files into a single FASTA file per sample.
    4. Executing the CAP3 assembly process for each sample.

    Args:
        id_map_file: Path to the ID map file (Excel format).
        raw_data_dir: Directory containing the raw sequencing files (e.g., .seq, .ab1).
        split_data_dir: Directory to store intermediate, sample-specific data.
        assembly_dir: Directory where the final assembly outputs will be stored.
    """
    logger.info("Loading ID map...")
    id_map_dict = data_handler.read_id_map(id_map_file)
    sample_list = list(id_map_dict.values())
    logger.info(f"Found {len(sample_list)} samples.")

    logger.info("Splitting raw data...")
    data_handler.split_data(id_map_dict, raw_data_dir, split_data_dir)

    logger.info("Combining sequences into FASTA files...")
    data_handler.combine_seq2fasta(split_data_dir, assembly_dir)

    _run_cap3_for_samples(assembly_dir, sample_list)


def main():
    """
    Parse command-line arguments and execute the genome assembly workflow.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="VSAT: Assemble raw sequencing data using CAP3.",
    )

    # 1. Input Data
    input_group = parser.add_argument_group("Input Data")
    input_group.add_argument(
        "-i", "--id_map", required=True, help="Path to the sample ID map file."
    )
    input_group.add_argument(
        "-r",
        "--raw_data",
        required=True,
        help="Directory containing raw sequencing files (.seq, .ab1).",
    )

    # 2. Workspace Options
    workspace_group = parser.add_argument_group("Workspace Options")
    workspace_group.add_argument(
        "-s",
        "--split_data",
        required=True,
        help="Directory to store intermediate organized split data.",
    )
    workspace_group.add_argument(
        "-a",
        "--assembly_dir",
        required=True,
        help="Directory where assembled contigs will be stored.",
    )

    args = parser.parse_args()

    run_genome_assembly(
        id_map_file=args.id_map,
        raw_data_dir=args.raw_data,
        split_data_dir=args.split_data,
        assembly_dir=args.assembly_dir,
    )


if __name__ == "__main__":
    main()

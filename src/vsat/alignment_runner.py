#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
This module provides a wrapper for running the MAFFT multiple sequence alignment tool.

It automates the process of preparing input files, executing the MAFFT command,
and then colorizing the resulting alignment for better visualization by calling
the `msa_colorizer.py` script.
"""

import sys
import subprocess
from pathlib import Path

from . import data_handler


def run_mafft_alignment(
    seq_dir: Path, out_dir: Path, alignment_type: str, ref_name: str
) -> Path | None:
    """
    Run a MAFFT alignment workflow for a set of sequences in a directory.

    This function performs the following steps:
    1. Merges all sequence files (.fasta, .fa, .seq) from `seq_dir` into a
       single multi-FASTA file, ensuring the reference sequence is listed first.
    2. Runs the MAFFT alignment tool on the merged file.
    3. Calls the `msa_colorizer.py` script to convert the alignment output
       into a color-coded HTML file.

    Note:
        This function requires the 'mafft' command-line tool to be installed
        and accessible in the system's PATH.

    Args:
        seq_dir: The directory containing the input sequence files.
        out_dir: The directory where all output files (merged input, alignment,
                 and colorized HTML) will be saved.
        alignment_type: A string used to name the output files (e.g., "genome").
        ref_name: The filename stem of the reference sequence, used to ensure
                  it appears first in the alignment.

    Returns:
        The path to the final colorized HTML alignment file if successful,
        otherwise None.
    """
    data_handler.mkdir(out_dir)

    merged_input_file = out_dir / f"{alignment_type}_merged.fasta"

    # Gather all supported sequence files from the input directory.
    all_seq_files = (
        list(seq_dir.glob("*.fa")) +
        list(seq_dir.glob("*.fasta")) +
        list(seq_dir.glob("*.seq"))
    )

    # Separate the reference sequence from the others to control its order.
    ref_seq_files = [f for f in all_seq_files if f.stem == ref_name]
    other_seq_files = [f for f in all_seq_files if f.stem != ref_name]

    # Recombine the file list, placing the reference first, followed by sorted others.
    seq_files = ref_seq_files + sorted(other_seq_files, key=lambda p: p.name)
    # print(f"seq_files: {seq_files}")

    if not seq_files:
        print(f"No sequence files found in {seq_dir} for alignment.")
        return None

    # Merge all sequences into a single FASTA file for MAFFT input.
    with merged_input_file.open("w", encoding="utf-8") as outfile:
        for seq_file in seq_files:
            seq_id, sequence = data_handler.read_genome_sequence(seq_file)
            outfile.write(f">{seq_id}\n{sequence}\n")

    # Run MAFFT alignment.
    mafft_output_file = out_dir / f"{alignment_type}_aligned.fasta"
    mafft_cmd = ["mafft", "--auto", str(merged_input_file)]

    print(f"Running MAFFT for {alignment_type}... Command: {' '.join(mafft_cmd)}")
    try:
        with mafft_output_file.open("w", encoding="utf-8") as f_out:
            result = subprocess.run(
                mafft_cmd,
                capture_output=True,
                text=True,
                check=True,
                encoding="utf-8",
            )
            f_out.write(result.stdout)
    except FileNotFoundError:
        print("Error: 'mafft' command not found. Please ensure it is installed and in your PATH.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error during MAFFT execution for {alignment_type}:")
        print(e.stderr)
        return None

    # Colorize the MAFFT output using the msa_colorizer.py script.
    colorized_output_html = out_dir / f"{alignment_type}_aligned_colorized.html"
    script_path = Path(__file__).resolve().parent / "msa_colorizer.py"
    if not script_path.exists():
        print(f"Error: Colorizing script not found at {script_path}")
        return None

    colorize_cmd = [
        sys.executable,
        str(script_path),
        str(mafft_output_file),
        "-o",
        str(colorized_output_html),
    ]

    print(f"Colorizing MAFFT output for {alignment_type}...")
    try:
        subprocess.run(
            colorize_cmd,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
    except subprocess.CalledProcessError as e:
        print(f"Error during colorizing script execution for {alignment_type}:")
        print(e.stderr)
        return None

    print(f"Successfully generated colorized alignment: {colorized_output_html}")
    return colorized_output_html

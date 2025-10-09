#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import subprocess
import sys
from pathlib import Path


def run_mafft_alignment(
    seq_dir: Path, out_dir: Path, alignment_type: str
) -> Path | None:
    """
    Runs MAFFT alignment on sequences found in a directory, merges them,
    and colorizes the output.
    """
    out_dir.mkdir(exist_ok=True)

    merged_input_file = out_dir / f"{alignment_type}_merged.fasta"
    seq_files = list(seq_dir.glob("*.fa*")) + list(seq_dir.glob("*.seq"))

    if not seq_files:
        print(f"No sequence files found in {seq_dir} for alignment.")
        return None

    with merged_input_file.open("w", encoding="utf-8") as outfile:
        for seq_file in seq_files:
            outfile.write(seq_file.read_text(encoding="utf-8"))

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
        print("Error: 'mafft' command not found. Please ensure it's installed and in your PATH.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error during MAFFT execution for {alignment_type}:")
        print(e.stderr)
        return None

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

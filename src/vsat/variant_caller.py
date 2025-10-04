#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
This script serves as the main entry point for the SNP detection workflow.
It coordinates the reading of data, alignment, and report generation.
"""
import argparse
from pathlib import Path
from . import data_handler, snp_analyzer, report_generator


# Define paths relative to the script location for default CLI behavior
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent.parent / "data"
GENOME_DIR = DATA_DIR / "genomes"
REF_GENOME_FILE = DATA_DIR / "ref_genome.json"


def _print_reference_genomes(genome_dict: dict[str, list[str]]) -> None:
    """Prints the available reference genomes in a formatted way."""
    print("Available reference genomes:")
    for virus_name, loci in genome_dict.items():
        print(f"{virus_name:<10}{', '.join(loci)}")


def run_snp_analysis(
    ref_genome_fasta_file: str | Path,
    ref_genome_gff_file: str | Path,
    assembled_sequences_dir: str | Path
) -> None:
    """
    Orchestrates the entire SNP analysis workflow.

    Args:
        assembled_sequences_dir: Directory containing assembled sample sequences.
        ref_genome_fasta_file: Path to the reference genome FASTA file.
        ref_genome_gff_file: Path to the reference genome GFF3 annotation file.
    """
    print("Starting SNP analysis...")

    print(f"Loading reference genome from: {ref_genome_fasta_file}")
    ref_name, ref_seq = data_handler.read_genome_sequence(ref_genome_fasta_file)

    print(f"Loading reference genome annotation from: {ref_genome_gff_file}")
    ref_peptides_dict = data_handler.parse_gff(ref_genome_gff_file)

    print(f"Loading assembled sequences from: {assembled_sequences_dir}")
    genome_seq_dict = data_handler.read_assembled_sequences(assembled_sequences_dir)
    print(f"Found {len(genome_seq_dict)} samples to analyze.")

    print("Aligning sequences and identifying mutations...")
    aligned_genome_seq_dict, mutate_sites_list = (
        snp_analyzer.align_sequences_and_find_mutations(ref_seq, genome_seq_dict)
    )

    print("Generating reports...")
    report_generator.write_snp_results_to_txt(
        aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq
    )
    print(" - snp_result.txt created.")
    report_generator.write_snp_results_to_xls(
        aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq
    )
    print(" - snp_result.xls created.")
    report_generator.align_stat(
        ref_peptides_dict, aligned_genome_seq_dict, ref_name, ref_seq
    )
    print(" - 测序比对结果统计.xls created.")

    print("Storing processed sequences...")
    # Add the reference sequence to the dictionary for inclusion in the output
    aligned_genome_seq_dict[ref_name] = ref_seq
    report_generator.store_sequence(
        assembled_sequences_dir,
        aligned_genome_seq_dict,
        ref_peptides_dict,
        ref_genome_fasta_file,
    )
    print(" - Processed sequences stored in '序列拼接' directory.")

    print("Analysis complete.")


def main():
    """
    Parses command-line arguments and executes the SNP detection workflow.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Analyzes viral SNPs from assembled sequences.",
    )
    parser.add_argument(
        "-g",
        "--locus",
        help=f"The locus of the virus to analyze. See available locus in: {GENOME_DIR}",
    )
    parser.add_argument(
        "-s",
        "--assembly_dir",
        help="Directory containing the assembled sequences (as .seq files).",
    )
    parser.add_argument(
        "-l",
        "--list",
        action="store_true",
        help="List available reference genomes and exit.",
    )

    args = parser.parse_args()

    if args.list:
        genome_dict = data_handler.load_reference_genomes(REF_GENOME_FILE)
        _print_reference_genomes(genome_dict)
        exit()

    if not args.locus or not args.assembly_dir:
        parser.error(
            "Both --locus and --assembly_dir are required unless --list is specified."
        )

    locus = args.locus.strip()
    assembly_dir = Path(args.assembly_dir)
    ref_genome_file = GENOME_DIR / f"{locus}.fasta"
    ref_genome_ann_file = GENOME_DIR / f"{locus}.gff3"

    if not ref_genome_file.is_file():
        parser.error(f"Genome file not found: {ref_genome_file}")
    if not ref_genome_ann_file.is_file():
        parser.error(f"Annotation file not found: {ref_genome_ann_file}")

    run_snp_analysis(
        ref_genome_fasta_file=ref_genome_file,
        ref_genome_gff_file=ref_genome_ann_file,
        assembled_sequences_dir=assembly_dir
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
This script serves as the main entry point for the SNP detection workflow.
It coordinates the reading of data, alignment, and report generation.
"""
import argparse
import os
from pathlib import Path
from . import data_handler, snp_analyzer, report_generator, alignment_runner


# Define paths relative to the script location for default CLI behavior
SCRIPT_DIR: Path = Path(__file__).resolve().parent
DATA_DIR: Path = SCRIPT_DIR.parent.parent / "data"
GENOME_DIR: Path = DATA_DIR / "genomes"
REF_GENOME_FILE: Path = DATA_DIR / "ref_genome.json"


def _print_reference_genomes(genome_dict: dict[str, list[str]]) -> None:
    """Prints the available reference genomes in a formatted way."""
    print("Available reference genomes:")
    for virus_name, loci in genome_dict.items():
        print(f"{virus_name:<10}{', '.join(loci)}")


def _phase_1_load_and_analyze(
    ref_genome_fasta_file: Path,
    ref_genome_gff_file: Path,
    assembled_sequences_dir: Path,
):
    """Loads all data and performs the core SNP analysis."""
    print("Step 1: Loading data and performing core analysis...")
    print(f" - Loading reference genome from: {ref_genome_fasta_file}")
    ref_name, ref_seq = data_handler.read_genome_sequence(ref_genome_fasta_file)
    print(f" - Loading reference genome annotation from: {ref_genome_gff_file}")
    ref_peptides_dict = data_handler.parse_gff(ref_genome_gff_file)
    print(f" - Loading assembled sequences from: {assembled_sequences_dir}")
    genome_seq_dict = data_handler.read_assembled_sequences(assembled_sequences_dir)
    print(f" - Found {len(genome_seq_dict)} samples to analyze.")

    print(" - Aligning sequences and identifying mutations...")
    aligned_genome_seq_dict, mutate_sites_list = (
        snp_analyzer.align_sequences_and_find_mutations(ref_seq, genome_seq_dict)
    )
    print("Step 1 complete.")
    return (
        ref_name,
        ref_seq,
        ref_peptides_dict,
        aligned_genome_seq_dict,
        mutate_sites_list,
    )


def _phase_2_generate_intermediate_reports(
    aligned_genome_seq_dict: dict[str, str],
    mutate_sites_list: list[str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_seq: str,
    ref_name: str,
):
    """Generates intermediate text and XLS reports."""
    print("\nStep 2: Generating intermediate text and XLS reports...")
    report_generator.align_stat(
        ref_peptides_dict, aligned_genome_seq_dict, ref_name, ref_seq
    )
    print(" - alignment_statistics.xls created.")

    report_generator.write_snp_results_to_txt(
        aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq
    )
    print(" - snp_result.txt created.")

    report_generator.write_snp_results_to_xls(
        aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq
    )
    print(" - snp_result.xls created.")
    print("Step 2 complete.")


def _phase_3_store_sequences(
    assembled_sequences_dir: Path,
    aligned_genome_seq_dict: dict[str, str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_genome_fasta_file: Path,
    ref_name: str,
    ref_seq: str,
    out_seq_dir: Path,
):
    """Stores processed sequences for the alignment phase."""
    print("\nStep 3: Storing processed sequences...")
    # Add the reference sequence to the dictionary for inclusion in the output files
    aligned_genome_seq_dict_with_ref = aligned_genome_seq_dict.copy()
    aligned_genome_seq_dict_with_ref[ref_name] = ref_seq
    report_generator.store_sequence(
        assembled_sequences_dir,
        aligned_genome_seq_dict_with_ref,
        ref_peptides_dict,
        ref_genome_fasta_file,
        out_seq_dir,
    )
    print(f" - Sequences stored in {out_seq_dir} directory.")
    print("Step 3 complete.")


def _phase_4_run_alignments(out_seq_dir: Path):
    """Runs MAFFT alignments for genome, all genes, and all proteins."""
    print("\nStep 4: Running MAFFT alignments for all sequences...")

    # Genome alignment
    genome_dir = out_seq_dir / "genome"
    if genome_dir.is_dir():
        alignment_runner.run_mafft_alignment(
            seq_dir=genome_dir,
            out_dir=genome_dir,  # Save alignment in the same directory
            alignment_type="genome",
        )

    # Gene alignments
    gene_seq_parent_dir = out_seq_dir / "gene"
    if gene_seq_parent_dir.is_dir():
        for subdir in sorted(gene_seq_parent_dir.iterdir()):
            if subdir.is_dir():
                alignment_runner.run_mafft_alignment(
                    seq_dir=subdir,
                    out_dir=subdir,  # Save alignment in the same directory
                    alignment_type=f"gene_{subdir.name}",
                )

    # Protein alignments
    protein_seq_parent_dir = out_seq_dir / "protein"
    if protein_seq_parent_dir.is_dir():
        for subdir in sorted(protein_seq_parent_dir.iterdir()):
            if subdir.is_dir():
                alignment_runner.run_mafft_alignment(
                    seq_dir=subdir,
                    out_dir=subdir,  # Save alignment in the same directory
                    alignment_type=f"protein_{subdir.name}",
                )
    print("Step 4 complete.")


def _phase_5_generate_final_report(out_seq_dir: Path):
    """Selects featured alignments and generates the final HTML report."""
    print("\nStep 5: Generating final HTML report...")
    # Select the first gene and protein alignment to feature in the report
    target_gene_dir = next(
        (d for d in (out_seq_dir / "gene").iterdir() if d.name.startswith("1_")),
        None,
    )
    target_protein_dir = next(
        (
            d
            for d in (out_seq_dir / "protein").iterdir()
            if d.name.startswith("1_")
        ),
        None,
    )

    # Construct paths to the alignment files in their new locations
    gene_mafft_file = (
        target_gene_dir / f"gene_{target_gene_dir.name}_aligned_colorized.html"
        if target_gene_dir
        else None
    )
    protein_mafft_file = (
        target_protein_dir / f"protein_{target_protein_dir.name}_aligned_colorized.html"
        if target_protein_dir
        else None
    )

    report_generator.generate_html_report(
        gene_mafft_file=gene_mafft_file,
        protein_mafft_file=protein_mafft_file,
        gene_dir_name=target_gene_dir.name if target_gene_dir else None,
        protein_dir_name=target_protein_dir.name if target_protein_dir else None,
    )
    print(" - viral_snp_report.html created.")
    print("Step 5 complete.")


def run_snp_analysis(
    ref_genome_fasta_file: Path,
    ref_genome_gff_file: Path,
    assembled_sequences_dir: Path,
    output_dir: Path,
) -> None:
    """
    Orchestrates the entire SNP analysis workflow by executing five distinct phases.

    Args:
        ref_genome_fasta_file: Path to the reference genome FASTA file.
        ref_genome_gff_file: Path to the reference genome GFF3 annotation file.
        assembled_sequences_dir: Directory containing assembled sample sequences.
        output_dir: The main directory to store all results.
    """
    print("Starting SNP analysis workflow...")
    out_seq_dir = Path("processed_sequences")

    # Phase 1: Data Loading and Core Analysis
    (
        ref_name,
        ref_seq,
        ref_peptides_dict,
        aligned_genome_seq_dict,
        mutate_sites_list,
    ) = _phase_1_load_and_analyze(
        ref_genome_fasta_file, ref_genome_gff_file, assembled_sequences_dir
    )

    # Phase 2: Generate Intermediate Reports
    _phase_2_generate_intermediate_reports(
        aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq, ref_name
    )

    # Phase 3: Store Processed Sequences for Alignment
    _phase_3_store_sequences(
        assembled_sequences_dir,
        aligned_genome_seq_dict,
        ref_peptides_dict,
        ref_genome_fasta_file,
        ref_name,
        ref_seq,
        out_seq_dir=out_seq_dir,
    )

    # Phase 4: Run MAFFT Alignments
    _phase_4_run_alignments(out_seq_dir=out_seq_dir)

    # Phase 5: Generate Final HTML Report
    _phase_5_generate_final_report(out_seq_dir=out_seq_dir)

    print("\nAnalysis complete.")


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
        "-o",
        "--output_dir",
        default="viral_snp_results",
        help="Directory to store all the results.",
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
    assembly_dir = Path(args.assembly_dir).resolve()
    output_dir = Path(args.output_dir)

    # Create the output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Change the current working directory to the output directory
    os.chdir(output_dir)
    print(f"Changed working directory to: {output_dir}")

    ref_genome_file = GENOME_DIR / f"{locus}.fasta"
    ref_genome_ann_file = GENOME_DIR / f"{locus}.gff3"

    if not ref_genome_file.is_file():
        parser.error(f"Genome file not found: {ref_genome_file}")
    if not ref_genome_ann_file.is_file():
        parser.error(f"Annotation file not found: {ref_genome_ann_file}")

    run_snp_analysis(
        ref_genome_fasta_file=ref_genome_file,
        ref_genome_gff_file=ref_genome_ann_file,
        assembled_sequences_dir=assembly_dir,
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()

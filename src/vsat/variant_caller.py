#!/usr/bin/env python3
# -*- coding:utf-8 -*-
"""
This script is the main entry point for the VSAT SNP detection workflow.

It orchestrates the entire analysis pipeline, which is divided into five phases:
1. Data Loading and Core Analysis: Loads reference and sample sequences,
   performs alignment, and identifies mutations.
2. Intermediate Report Generation: Creates text and Excel-compatible reports
   for alignment statistics and SNP details.
3. Sequence Storage: Organizes and saves processed sequences (full genome,
   genes, and proteins) into a structured directory for subsequent alignment.
4. MAFFT Alignment: Runs multiple sequence alignments for all sequence sets.
5. Final HTML Report Generation: Compiles all results into a single,
   comprehensive HTML report with embedded visualizations.
"""

import argparse
import logging
import sys
import shutil
from pathlib import Path

from concurrent.futures import ThreadPoolExecutor
from . import (
    data_handler,
    snp_analyzer,
    report_generator,
    alignment_runner,
    logger as vsat_logger,
)


# Set up logging with colored output.
logger = vsat_logger.setup_logger(__name__)


# Define default paths relative to the script's location.
SCRIPT_DIR: Path = Path(__file__).resolve().parent
DATA_DIR: Path = SCRIPT_DIR.parent.parent / "data"
GENOME_DIR: Path = DATA_DIR / "genomes"
REF_GENOME_FILE: Path = DATA_DIR / "ref_genome.json"


def _print_reference_genomes(genome_dict: dict[str, list[str]]) -> None:
    """Print the available reference genomes in a formatted list."""
    print("Available reference genomes:")
    for virus_name, loci in genome_dict.items():
        print(f"{virus_name:<10}{', '.join(loci)}")


def _phase_1_load_and_analyze(
    ref_genome_seq_file: Path,
    ref_genome_gff_file: Path,
    assembled_sequences_dir: Path,
) -> tuple[str, str, dict[str, tuple[int, int]], dict[str, str], list[str]]:
    """Load all necessary data and perform the core SNP analysis."""
    logger.info("Step 1: Loading data and performing core analysis...")
    logger.info(f" - Reference genome:\n   {ref_genome_seq_file}")
    ref_name, ref_seq = data_handler.read_genome_sequence(seq_file=ref_genome_seq_file)
    logger.info(f" - Reference annotation:\n   {ref_genome_gff_file}")
    ref_peptides_dict = data_handler.parse_gff(gff3_file=ref_genome_gff_file)
    logger.info(f" - Assembled sequences:\n   {assembled_sequences_dir}")
    genome_seq_dict = data_handler.read_assembled_sequences(
        directory=assembled_sequences_dir
    )
    logger.info(f" - Found {len(genome_seq_dict)} samples to analyze.")

    logger.info(" - Aligning sequences and identifying mutations...")
    aligned_genome_seq_dict, mutate_sites_list = (
        snp_analyzer.align_sequences_and_find_mutations(
            ref_seq=ref_seq, genome_seq_dict=genome_seq_dict
        )
    )
    logger.info("Step 1 complete.")
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
    output_dir: Path,
):
    """Generate intermediate analysis reports in text and XLS formats."""
    logger.info("\nStep 2: Generating intermediate text and XLS reports...")
    report_generator.align_stat(
        ref_peptides_dict=ref_peptides_dict,
        aligned_genome_seq_dict=aligned_genome_seq_dict,
        ref_name=ref_name,
        ref_seq=ref_seq,
        output_dir=output_dir,
    )
    logger.info(f" - alignment_statistics.xls created in:\n   {output_dir}")

    report_generator.write_snp_results_to_txt(
        aligned_genome_seq_dict=aligned_genome_seq_dict,
        mutate_sites_list=mutate_sites_list,
        ref_peptides_dict=ref_peptides_dict,
        ref_seq=ref_seq,
        output_dir=output_dir,
    )
    logger.info(f" - snp_result.txt created in:\n   {output_dir}")

    report_generator.write_snp_results_to_xls(
        aligned_genome_seq_dict=aligned_genome_seq_dict,
        mutate_sites_list=mutate_sites_list,
        ref_peptides_dict=ref_peptides_dict,
        ref_seq=ref_seq,
        output_dir=output_dir,
    )
    logger.info(f" - snp_result.xls created in:\n   {output_dir}")
    logger.info("Step 2 complete.")


def _phase_3_store_sequences(
    assembled_sequences_dir: Path,
    aligned_genome_seq_dict: dict[str, str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_genome_seq_file: Path,
    ref_name: str,
    ref_seq: str,
    out_seq_dir: Path,
):
    """Store processed sequences in a structured directory for the alignment phase."""
    logger.info("\nStep 3: Storing processed sequences...")
    # Include the reference sequence in the dictionary to ensure it's part of the output.
    aligned_genome_seq_dict_with_ref = aligned_genome_seq_dict.copy()
    aligned_genome_seq_dict_with_ref[ref_name] = ref_seq
    report_generator.store_sequence(
        assembly_dir=assembled_sequences_dir,
        aligned_genome_seq_dict=aligned_genome_seq_dict_with_ref,
        ref_peptides_dict=ref_peptides_dict,
        ref_genome_file=ref_genome_seq_file,
        out_seq_dir=out_seq_dir,
    )
    logger.info(f" - Sequences stored in:\n   {out_seq_dir}")
    logger.info("Step 3 complete.")


def _phase_4_run_alignments(out_seq_dir: Path, ref_name: str, max_workers: int = 4):
    """Run MAFFT alignments for the full genome, all genes, and all proteins in parallel."""
    logger.info("\nStep 4: Running MAFFT alignments for all sequences in parallel...")

    tasks = []
    # Run genome alignment
    genome_dir = out_seq_dir / "genome"
    if genome_dir.is_dir():
        tasks.append((genome_dir, genome_dir, "genome", ref_name))

    # Run gene alignments for each gene subdirectory
    gene_seq_parent_dir = out_seq_dir / "gene"
    if gene_seq_parent_dir.is_dir():
        for subdir in sorted(gene_seq_parent_dir.iterdir()):
            if subdir.is_dir():
                tasks.append((subdir, subdir, f"gene_{subdir.name}", ref_name))

    # Run protein alignments for each protein subdirectory
    protein_seq_parent_dir = out_seq_dir / "protein"
    if protein_seq_parent_dir.is_dir():
        for subdir in sorted(protein_seq_parent_dir.iterdir()):
            if subdir.is_dir():
                tasks.append((subdir, subdir, f"protein_{subdir.name}", ref_name))

    # Use ThreadPoolExecutor to run alignments in parallel.
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(alignment_runner.run_mafft_alignment, *task)
            for task in tasks
        ]
        # Wait for all tasks to complete.
        for future in futures:
            future.result()

    logger.info("Step 4 complete.")


def _phase_5_generate_final_report(out_seq_dir: Path, output_dir: Path):
    """Select featured alignments and generate the final, comprehensive HTML report."""
    logger.info("\nStep 5: Generating final HTML report...")
    # Select the first gene and protein alignment (sorted alphabetically) to feature.
    target_gene_dir = next(
        (
            d
            for d in (out_seq_dir / "gene").iterdir()
            if d.is_dir() and d.name.startswith("1.")
        ),
        None,
    )
    target_protein_dir = next(
        (
            d
            for d in (out_seq_dir / "protein").iterdir()
            if d.is_dir() and d.name.startswith("1.")
        ),
        None,
    )

    # Construct paths to the colorized alignment files to be embedded.
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
        output_dir=output_dir,
    )
    logger.info(f" - viral_snp_report.html created in:\n   {output_dir}")
    logger.info("Step 5 complete.")


def run_snp_analysis(
    ref_genome_seq_file: Path,
    ref_genome_gff_file: Path,
    assembled_sequences_dir: Path,
    output_dir: Path,
    max_workers: int = 4,
) -> None:
    """
    Orchestrate the entire SNP analysis workflow through its five phases.

    Args:
        ref_genome_seq_file: Path to the reference genome FASTA file.
        ref_genome_gff_file: Path to the reference genome GFF3 annotation file.
        assembled_sequences_dir: Directory containing assembled sample sequences.
        output_dir: The main directory where all results will be stored.
        max_workers: The maximum number of threads to use for parallel tasks.
    """
    logger.info("Starting SNP analysis workflow...")
    out_seq_dir = output_dir / "processed_sequences"

    # Phase 1: Data Loading and Core Analysis
    (
        ref_name,
        ref_seq,
        ref_peptides_dict,
        aligned_genome_seq_dict,
        mutate_sites_list,
    ) = _phase_1_load_and_analyze(
        ref_genome_seq_file=ref_genome_seq_file,
        ref_genome_gff_file=ref_genome_gff_file,
        assembled_sequences_dir=assembled_sequences_dir,
    )

    # Phase 2: Generate Intermediate Reports
    _phase_2_generate_intermediate_reports(
        aligned_genome_seq_dict=aligned_genome_seq_dict,
        mutate_sites_list=mutate_sites_list,
        ref_peptides_dict=ref_peptides_dict,
        ref_seq=ref_seq,
        ref_name=ref_name,
        output_dir=output_dir,
    )

    # Phase 3: Store Processed Sequences for Alignment
    _phase_3_store_sequences(
        assembled_sequences_dir=assembled_sequences_dir,
        aligned_genome_seq_dict=aligned_genome_seq_dict,
        ref_peptides_dict=ref_peptides_dict,
        ref_genome_seq_file=ref_genome_seq_file,
        ref_name=ref_name,
        ref_seq=ref_seq,
        out_seq_dir=out_seq_dir,
    )

    # Phase 4: Run MAFFT Alignments in parallel
    _phase_4_run_alignments(
        out_seq_dir=out_seq_dir, ref_name=ref_name, max_workers=max_workers
    )

    # Phase 5: Generate Final HTML Report
    _phase_5_generate_final_report(out_seq_dir=out_seq_dir, output_dir=output_dir)

    logger.info("\nAnalysis complete.")


def main():
    """
    Parse command-line arguments and execute the SNP detection workflow.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="VSAT: Analyzes viral SNPs from assembled sequences.",
    )

    # 1. Core Inputs
    input_group = parser.add_argument_group("Core Inputs")
    input_group.add_argument(
        "-g",
        "--locus",
        help=f"The locus of the virus to analyze. See available locus in: {GENOME_DIR}",
    )
    input_group.add_argument(
        "-s",
        "--assembly_dir",
        help="Directory containing the assembled sequences (as .seq files).",
    )

    # 2. Output Options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument(
        "-o",
        "--output_dir",
        default="viral_snp_results",
        help="Directory to store all the results.",
    )

    # 3. Performance Options
    perf_group = parser.add_argument_group("Performance Options")
    perf_group.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Maximum number of threads for parallel alignments.",
    )

    # 4. Utility
    utility_group = parser.add_argument_group("Utility")
    utility_group.add_argument(
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
    output_dir = Path(args.output_dir).resolve()

    # Clear the output directory if it already exists to ensure a fresh start.
    if output_dir.exists():
        logger.info(f"Clearing existing output directory: {output_dir}")
        shutil.rmtree(output_dir)

    # Create the output directory.
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_genome_file = GENOME_DIR / f"{locus}.fasta"
    ref_genome_ann_file = GENOME_DIR / f"{locus}.gff3"

    if not ref_genome_file.is_file():
        parser.error(f"Genome file not found: {ref_genome_file}")
    if not ref_genome_ann_file.is_file():
        parser.error(f"Annotation file not found: {ref_genome_ann_file}")

    run_snp_analysis(
        ref_genome_seq_file=ref_genome_file,
        ref_genome_gff_file=ref_genome_ann_file,
        assembled_sequences_dir=assembly_dir,
        output_dir=output_dir,
        max_workers=args.threads,
    )


if __name__ == "__main__":
    main()

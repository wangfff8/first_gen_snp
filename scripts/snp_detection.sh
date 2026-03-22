#!/bin/bash

# 2. SNP Detection
# Description: This script runs the SNP (Single Nucleotide Polymorphism) detection workflow.
# It first lists the available reference genomes and then provides a template command
# to run the actual SNP analysis using the `vsat.variant_caller` module.
#
# --- Step 1: List Available Reference Genomes ---
# This command shows the shorthand locus names for the reference genomes
# defined in `data/ref_genome.json`. You will need one of these names for the `-g` parameter below.
echo "INFO: Listing available reference genomes..."
python -m vsat.variant_caller --list
echo "----------------------------------------"
echo ""

# --- Step 2: Run SNP Detection ---
# The command to run the SNP detection is provided below.
#
# IMPORTANT:
# 1. Uncomment the `python -m vsat.variant_caller ...` line below.
# 2. Replace `locus_name` with a valid locus name from the list above.
# 3. Replace `/path/to/final_assembly_dir` with the path to your assembled sequences.
# 4. Replace `/path/to/output_results_dir` with the path where you want to save the results.
# 5. Use `-t` or `--threads` to specify the number of threads for parallel alignments.
#
echo "INFO: To run SNP detection, please edit this script and uncomment the following line."
echo "INFO: Replace the placeholder values with your specific data."

echo "# python -m vsat.variant_caller -g <locus_name> -s /path/to/final_assembly_dir -o /path/to/output_results_dir -t 4"

echo "Script execution finished."

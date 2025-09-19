#!/bin/bash
set -e

# 1. Genome Assembly
# This script provides an example of how to run the genome assembly workflow.
# Before executing, ensure you have activated the virtual environment.
# For example, on Linux/macOS: source .venv/bin/activate
# For example, on Windows:   source .venv/Scripts/activate

# Replace the placeholder paths below with the actual paths to your data.
python -m vsat.genome_assembler \
    --rawdata /path/to/rawdata \
    --split_data /path/to/split_data \
    --assemble_dir /path/to/assemble \
    --id_map /path/to/id_map.xls

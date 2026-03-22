"""
VSAT: Virus SNP Analysis Tool

A Python package for analyzing Single Nucleotide Polymorphisms (SNPs) in viral
genomes using first-generation sequencing data.
"""

__version__ = "0.1.1"

from . import (
    data_handler,
    snp_analyzer,
    report_generator,
    alignment_runner,
)

__all__ = [
    "data_handler",
    "snp_analyzer",
    "report_generator",
    "alignment_runner",
]

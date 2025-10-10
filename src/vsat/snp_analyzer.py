#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
This module provides the core bioinformatics functions for SNP analysis.

It includes functionalities such as sequence alignment and padding,
mutation identification, sequence identity calculation, codon-to-amino-acid
translation, and detailed SNP effect analysis. The module relies on a
predefined codon table for translations.
"""

# A comprehensive mapping from DNA codons to amino acids.
# Each codon maps to a list containing:
#   - Single-letter amino acid code (e.g., 'A' for Alanine)
#   - Three-letter amino acid code (e.g., 'Ala')
#   - Full name of the amino acid in Chinese (e.g., '丙氨酸')
# Stop codons are represented with an asterisk ('*').
aa_dict: dict[str, list[str]] = {
    'GCT': ['A', 'Ala', '丙氨酸'], 'GCC': ['A', 'Ala', '丙氨酸'], 'GCA': ['A', 'Ala', '丙氨酸'],
    'GCG': ['A', 'Ala', '丙氨酸'], 'CGT': ['R', 'Arg', '精氨酸'], 'CGC': ['R', 'Arg', '精氨酸'],
    'CGA': ['R', 'Arg', '精氨酸'], 'CGG': ['R', 'Arg', '精氨酸'], 'AGA': ['R', 'Arg', '精氨酸'],
    'AGG': ['R', 'Arg', '精氨酸'], 'AAT': ['N', 'Asn', '天冬酰胺'], 'AAC': ['N', 'Asn', '天冬酰胺'],
    'GAT': ['D', 'Asp', '天冬氨酸'], 'GAC': ['D', 'Asp', '天冬氨酸'], 'TGT': ['C', 'Cys', '半胱氨酸'],
    'TGC': ['C', 'Cys', '半胱氨酸'], 'CAA': ['Q', 'Gln', '谷氨酰胺'], 'CAG': ['Q', 'Gln', '谷氨酰胺'],
    'GAA': ['E', 'Glu', '谷氨酸'], 'GAG': ['E', 'Glu', '谷氨酸'], 'GGT': ['G', 'Gly', '甘氨酸'],
    'GGC': ['G', 'Gly', '甘氨酸'], 'GGA': ['G', 'Gly', '甘氨酸'], 'GGG': ['G', 'Gly', '甘氨酸'],
    'CAT': ['H', 'His', '组氨酸'], 'CAC': ['H', 'His', '组氨酸'], 'ATT': ['I', 'Ile', '异亮氨酸'],
    'ATC': ['I', 'Ile', '异亮氨酸'], 'ATA': ['I', 'Ile', '异亮氨酸'], 'TTA': ['L', 'Leu', '亮氨酸'],
    'TTG': ['L', 'Leu', '亮氨酸'], 'CTT': ['L', 'Leu', '亮氨酸'], 'CTC': ['L', 'Leu', '亮氨酸'],
    'CTA': ['L', 'Leu', '亮氨酸'], 'CTG': ['L', 'Leu', '亮氨酸'], 'AAA': ['K', 'Lys', '赖氨酸'],
    'AAG': ['K', 'Lys', '赖氨酸'], 'ATG': ['M', 'Met', '甲硫氨酸'], 'TTT': ['F', 'Phe', '苯丙氨酸'],
    'TTC': ['F', 'Phe', '苯丙氨酸'], 'CCT': ['P', 'Pro', '脯氨酸'], 'CCC': ['P', 'Pro', '脯氨酸'],
    'CCA': ['P', 'Pro', '脯氨酸'], 'CCG': ['P', 'Pro', '脯氨酸'], 'TCT': ['S', 'Ser', '丝氨酸'],
    'TCC': ['S', 'Ser', '丝氨酸'], 'TCA': ['S', 'Ser', '丝氨酸'], 'TCG': ['S', 'Ser', '丝氨酸'],
    'AGT': ['S', 'Ser', '丝氨酸'], 'AGC': ['S', 'Ser', '丝氨酸'], 'ACT': ['T', 'Thr', '苏氨酸'],
    'ACC': ['T', 'Thr', '苏氨酸'], 'ACA': ['T', 'Thr', '苏氨酸'], 'ACG': ['T', 'Thr', '苏氨酸'],
    'TGG': ['W', 'Trp', '色氨酸'], 'TAT': ['Y', 'Tyr', '酪氨酸'], 'TAC': ['Y', 'Tyr', '酪氨酸'],
    'GTT': ['V', 'Val', '缬氨酸'], 'GTC': ['V', 'Val', '缬氨酸'], 'GTA': ['V', 'Val', '缬氨酸'],
    'GTG': ['V', 'Val', '缬氨酸'], 'TAA': ['*', '*', '终止密码子'], 'TAG': ['*', '*', '终止密码子'],
    'TGA': ['*', '*', '终止密码子'],
}


def align_and_pad_sequence(ref_seq: str, query_seq: str, n: int = 10) -> str:
    """
    Align a query sequence to a reference sequence.

    This function finds the starting position of the query sequence within the
    reference by matching an initial segment of `n` characters. It then pads
    the query sequence with 'N' characters at the beginning and end to match
    the total length of the reference sequence.

    Args:
        ref_seq: The reference DNA sequence.
        query_seq: The query DNA sequence to be aligned and padded.
        n: The length of the initial segment of the query sequence used to
           find the alignment start position. Defaults to 10.

    Returns:
        The query sequence, padded with 'N's to match the reference length.

    Raises:
        ValueError: If the initial segment of the query sequence cannot be
                    found in the reference sequence.
    """
    start_pos = ref_seq.find(query_seq[:n])
    if start_pos == -1:
        raise ValueError(f"Initial segment '{query_seq[:n]}' not found in reference.")

    padding_left = "N" * start_pos
    padding_right = "N" * (len(ref_seq) - (start_pos + len(query_seq)))

    return padding_left + query_seq + padding_right


def find_mutations(ref_seq: str, query_seq: str) -> str:
    """
    Identify single nucleotide mutations between two aligned sequences.

    Compares a reference sequence and a query sequence of the same length,
    character by character, to find differences. Bases 'N' in the query
    sequence are ignored.

    Args:
        ref_seq: The reference DNA sequence.
        query_seq: The aligned and padded query DNA sequence.

    Returns:
        A comma-separated string listing all mutations, where each mutation
        is formatted as "position_ref>query" (e.g., "10_A>G,25_C>T").
        Returns an empty string if no mutations are found.
    """
    mutations = []
    for i, (ref_base, query_base) in enumerate(zip(ref_seq, query_seq)):
        if query_base.upper() != 'N' and ref_base != query_base:
            mutations.append(f"{i + 1}_{ref_base}>{query_base}")
    return ",".join(mutations)


def cal_identity(seq1: str, seq2: str) -> str:
    """
    Calculate the sequence identity percentage between two sequences.

    Args:
        seq1: The first DNA sequence.
        seq2: The second DNA sequence, which must be of the same length as seq1.

    Returns:
        A string representing the identity percentage, formatted to two
        decimal places (e.g., "99.85%").
    """
    matches = sum(1 for x, y in zip(seq1, seq2) if x == y)
    return f"{(matches / len(seq1)):.2%}"


def get_aa_seq(nuc_seq: str) -> str:
    """
    Translate a nucleotide sequence into a single-letter amino acid sequence.

    The translation is based on the standard codon table defined in `aa_dict`.
    Incomplete codons at the end of the sequence are ignored.

    Args:
        nuc_seq: The nucleotide (DNA) sequence to be translated.

    Returns:
        The resulting single-letter amino acid sequence. Unknown codons are
        translated as '*'.
    """
    aa_seq = []
    # Ensure that only full codons are processed by truncating the sequence.
    end_pos = len(nuc_seq) - (len(nuc_seq) % 3)
    for i in range(0, end_pos, 3):
        codon = nuc_seq[i:i+3]
        aa_seq.append(aa_dict.get(codon, ['*'])[0])
    return "".join(aa_seq)


def align_sequences_and_find_mutations(
    ref_seq: str, genome_seq_dict: dict[str, str]
) -> tuple[dict[str, str], list[str]]:
    """
    Align multiple sample sequences to a reference and identify mutations.

    This function iterates through a dictionary of sample sequences, aligns each
    to the reference sequence using `align_and_pad_sequence`, and then finds
    mutations for each aligned sequence using `find_mutations`.

    Args:
        ref_seq: The reference DNA sequence.
        genome_seq_dict: A dictionary where keys are sample names and values
                         are the corresponding DNA sequences.

    Returns:
        A tuple containing two elements:
        - A dictionary of aligned sequences, with the same keys as the input.
        - A list of comma-separated mutation strings, one for each sample.
    """
    aligned_genome_seq_dict = {}
    mutate_sites_list = []

    for sample_name, sample_sequence in genome_seq_dict.items():
        aligned_sequence = align_and_pad_sequence(ref_seq, sample_sequence)
        aligned_genome_seq_dict[sample_name] = aligned_sequence
        mutate_sites_list.append(find_mutations(ref_seq, aligned_sequence))

    return aligned_genome_seq_dict, mutate_sites_list


def get_snp_info(
    nuc_site: int, ref_seq: str, mutated_genome_seq: str, coding_range: tuple[int, int]
) -> tuple[int, str, str, str, str]:
    """
    Analyze the effect of a Single Nucleotide Polymorphism (SNP) on translation.

    Given a nucleotide position of a SNP, this function determines its impact
    on the corresponding codon and the resulting amino acid.

    Args:
        nuc_site: The 1-based position of the SNP in the full genome sequence.
        ref_seq: The complete reference DNA sequence.
        mutated_genome_seq: The mutated genome sequence, aligned to the reference.
        coding_range: A tuple (start, end) indicating the 1-based start and
                      end positions of the coding sequence (CDS) within the
                      genome.

    Returns:
        A tuple containing:
        - aa_n (int): The 1-based position of the affected amino acid.
        - ori_codon (str): The original codon from the reference sequence.
        - after_codon (str): The mutated codon from the query sequence.
        - ori_pep (str): The original amino acid (Chinese name).
        - after_pep (str): The mutated amino acid (Chinese name).
    """
    # Calculate the SNP position relative to the start of the coding sequence (1-based).
    cod_site = nuc_site - coding_range[0] + 1

    # Determine the affected amino acid position within the protein (1-based).
    aa_n = (cod_site - 1) // 3 + 1

    # Find the start of the affected codon in the full genome sequence (0-based).
    codon_start_in_cds = (aa_n - 1) * 3
    codon_start_in_genome = coding_range[0] + codon_start_in_cds - 1

    # Extract the original and mutated codons.
    ori_codon = ref_seq[codon_start_in_genome : codon_start_in_genome + 3]
    after_codon = mutated_genome_seq[codon_start_in_genome : codon_start_in_genome + 3]

    # Translate codons to amino acids, using their Chinese names.
    ori_pep = aa_dict.get(ori_codon, ['*', '*', '*'])[2]
    after_pep = aa_dict.get(after_codon, ['*', '*', '*'])[2]

    return (aa_n, ori_codon, after_codon, ori_pep, after_pep)

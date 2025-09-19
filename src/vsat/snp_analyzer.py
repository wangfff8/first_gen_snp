#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from typing import Any


# Codon to Amino Acid mapping
# Format: [Single Letter, Three Letter, Chinese Name]
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
    Aligns a query sequence to a reference by finding its start position
    and pads the missing ends with 'N' characters to match the reference length.

    Args:
        ref_seq: The reference sequence.
        query_seq: The query sequence to align.
        n: The length of the initial segment of the query sequence to use for finding
           the start position.

    Returns:
        The padded query sequence.

    Raises:
        ValueError: If the initial segment of the query sequence is not found.
    """
    start_pos = ref_seq.find(query_seq[:n])
    if start_pos == -1:
        raise ValueError(f"Initial segment '{query_seq[:n]}' not found in reference.")

    padding_left = "N" * start_pos
    padding_right = "N" * (len(ref_seq) - (start_pos + len(query_seq)))

    return padding_left + query_seq + padding_right

def find_mutations(ref_seq: str, query_seq: str) -> str:
    """
    Compares two sequences and finds all single nucleotide mutations.

    Args:
        ref_seq: The reference sequence.
        query_seq: The aligned and padded query sequence.

    Returns:
        A comma-separated string of mutations, e.g., "10_A>G,25_C>T".
    """
    mutations = []
    for i, (ref_base, query_base) in enumerate(zip(ref_seq, query_seq)):
        if query_base.upper() != 'N' and ref_base != query_base:
            mutations.append(f"{i + 1}_{ref_base}>{query_base}")
    return ",".join(mutations)

def cal_identity(seq1: str, seq2: str) -> str:
    """Calculates the identity percentage between two sequences."""
    if not seq1:
        return "0.00%"
    matches = sum(1 for x, y in zip(seq1, seq2) if x == y)
    return f"{(matches / len(seq1)):.2%}"

def get_aa_seq(nuc_seq: str) -> str:
    """Translates a nucleotide sequence into a single-letter amino acid sequence."""
    aa_seq = []
    # Ensure we only process full codons
    end_pos = len(nuc_seq) - (len(nuc_seq) % 3)
    for i in range(0, end_pos, 3):
        codon = nuc_seq[i:i+3]
        aa_seq.append(aa_dict.get(codon, ['*'])[0])
    return "".join(aa_seq)

def align_sequences_and_find_mutations(
    ref_seq: str, genome_seq_dict: dict[str, str]
) -> tuple[dict[str, str], list[str]]:
    """
    Aligns a dictionary of sample sequences to a reference and finds mutations for each.
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
    For a given SNP, determines its effect on the amino acid translation.

    Returns:
        A tuple containing: (amino_acid_position, original_codon,
        mutated_codon, original_amino_acid, mutated_amino_acid).
    """
    # Position within the coding sequence (1-based)
    cod_site = nuc_site - coding_range[0] + 1

    # Determine the affected amino acid position (1-based)
    aa_n = (cod_site - 1) // 3 + 1

    # Determine the start of the codon in the genome sequence
    codon_start_in_cds = (aa_n - 1) * 3
    codon_start_in_genome = coding_range[0] + codon_start_in_cds - 1

    # Get original and mutated codons
    ori_codon = ref_seq[codon_start_in_genome : codon_start_in_genome + 3]
    after_codon = mutated_genome_seq[codon_start_in_genome : codon_start_in_genome + 3]

    # Translate to amino acids (using the Chinese name)
    ori_pep = aa_dict.get(ori_codon, ['*', '*', '*'])[2]
    after_pep = aa_dict.get(after_codon, ['*', '*', '*'])[2]

    return (aa_n, ori_codon, after_codon, ori_pep, after_pep)

#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import shutil
from pathlib import Path
from typing import Any
from . import data_handler, snp_analyzer


# Report templates for text output
templates: dict[int, str] = {
    1: "第 {1} 个碱基 ({2}) 发生突变, 不影响 {3} 区第 {4} 个氨基酸的表达 ({7}); ",
    2: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 发生突变, 不影响 {3} 区第 {4} 个氨基酸的表达 ({7}); ",
    3: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 及第 {10} 碱基 ({11}) 发生突变, 不影响 {3} 区第 {4} 个氨基酸的表达 ({7}); ",
    4: "第 {1} 个碱基 ({2}) 发生突变, 造成 {3} 区第 {4} 个氨基酸由{7}变成{8}; ",
    5: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 发生突变, 造成 {3} 区第 {4} 个氨基酸由{7}变成{8}; ",
    6: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 及第 {10} 碱基 ({11}) 发生突变, 造成 {3} 区第 {4} 个氨基酸由{7}变成{8}; ",
    7: "第 {} 个碱基 ({}) 发生突变; ",
}


def write_snp_results_to_txt(
    aligned_genome_seq_dict: dict[str, str],
    mutate_sites_list: list[str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_seq: str,
) -> None:
    """Writes a detailed SNP analysis report in plain text format."""
    with open("snp_result.txt", "w", encoding="utf-8") as res:
        for (sample, mutated_genome_seq), mutate_sites in zip(
            aligned_genome_seq_dict.items(), mutate_sites_list
        ):
            res.write(f"{sample}:\n")

            if not mutate_sites:
                continue

            mutates_dict = {
                int(m.split("_")[0]): m.split("_")[1] for m in mutate_sites.split(",")
            }
            mutated_sites = set(mutates_dict.keys())
            in_cds_sites = set()

            for nuc_site in mutated_sites:
                for product_name, coding_range in ref_peptides_dict.items():
                    if coding_range[0] <= nuc_site <= coding_range[1]:
                        in_cds_sites.add(nuc_site)
                        snp_info = (
                            sample,
                            nuc_site,
                            mutates_dict[nuc_site],
                            product_name,
                        ) + snp_analyzer.get_snp_info(
                            nuc_site, ref_seq, mutated_genome_seq, coding_range
                        )

                        template_key = 4 if snp_info[-1] != snp_info[-2] else 1
                        res.write(templates[template_key].format(*snp_info))

            out_cds_sites = sorted(list(mutated_sites - in_cds_sites))
            for nuc_site in out_cds_sites:
                res.write(templates[7].format(nuc_site, mutates_dict[nuc_site]))
            res.write("\n")


def write_snp_results_to_xls(
    aligned_genome_seq_dict: dict[str, str],
    mutate_sites_list: list[str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_seq: str,
) -> None:
    """Writes a tab-separated SNP analysis report compatible with Excel."""
    with open("snp_result.xls", "w", encoding="gbk") as content:
        # Define and write the header
        header = [
            "Sample",
            "Nuc Site",
            "Mutation",
            "Product",
            "AA Pos",
            "Ref Codon",
            "Sample Codon",
            "Ref AA",
            "Sample AA",
        ]
        content.write("\t".join(header) + "\n")
        for (sample, mutated_genome_seq), mutate_sites in zip(
            aligned_genome_seq_dict.items(), mutate_sites_list
        ):
            if not mutate_sites:
                continue

            mutates_dict = {
                int(m.split("_")[0]): m.split("_")[1] for m in mutate_sites.split(",")
            }
            mutated_sites = set(mutates_dict.keys())

            for nuc_site in sorted(list(mutated_sites)):
                in_cds = False
                for product_name, coding_range in ref_peptides_dict.items():
                    if coding_range[0] <= nuc_site <= coding_range[1]:
                        in_cds = True
                        snp_info = (
                            sample,
                            nuc_site,
                            mutates_dict[nuc_site],
                            product_name,
                        ) + snp_analyzer.get_snp_info(
                            nuc_site, ref_seq, mutated_genome_seq, coding_range
                        )
                        content.write("\t".join(map(str, snp_info)) + "\n")
                        break
                if not in_cds:
                    content.write(f"{sample}\t{nuc_site}\t{mutates_dict[nuc_site]}\n")


def align_stat(
    ref_peptides_dict: dict[str, tuple[int, int]],
    aligned_genome_seq_dict: dict[str, str],
    ref_name: str,
    ref_seq: str,
) -> None:
    """Calculates and writes alignment identity statistics to an Excel-compatible file."""
    with open("测序比对结果统计.xls", "w", encoding="gbk") as stat:
        pep_features = [ref_name] + list(ref_peptides_dict.keys())
        stat.write("\t".join(h for h in pep_features) + "\n")

        for sample, sample_seq in aligned_genome_seq_dict.items():
            stat.write(f"{sample}\t{snp_analyzer.cal_identity(ref_seq, sample_seq)}\t")
            pep_align_rate = []
            for pep_id, coding_range in ref_peptides_dict.items():
                ref_nuc_seq = ref_seq[coding_range[0] - 1 : coding_range[1]]
                qry_nuc_seq = sample_seq[coding_range[0] - 1 : coding_range[1]]
                ref_aa_seq = snp_analyzer.get_aa_seq(ref_nuc_seq)
                qry_aa_seq = snp_analyzer.get_aa_seq(qry_nuc_seq)
                pep_align_rate.append(snp_analyzer.cal_identity(ref_aa_seq, qry_aa_seq))
            stat.write("\t".join(pep_align_rate) + "\n")


def store_sequence(
    assemble_dir: str | Path,
    aligned_genome_seq_dict: dict[str, str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_genome_file: str | Path,
    out_seq_dir: str | Path = "序列拼接",
) -> None:
    """
    Stores the processed sequences, categorized by product, into a directory structure.
    """
    assemble_dir = Path(assemble_dir)
    out_seq_dir = Path(out_seq_dir)
    ref_genome_file = Path(ref_genome_file)
    data_handler.mkdir(out_seq_dir)

    # Copy reference and sample genomes
    genome_dir = out_seq_dir / "genome"
    data_handler.mkdir(genome_dir)
    if ref_genome_file.exists():
        shutil.copy(ref_genome_file, genome_dir)

    for sample_seq_file in assemble_dir.glob("*.seq"):
        if sample_seq_file.exists():
            shutil.copy(sample_seq_file, genome_dir)

    # Create categorized nucleotide and protein sequence files
    nuc_dir = out_seq_dir / "nucleotide"
    pep_dir = out_seq_dir / "protein"
    data_handler.mkdir(nuc_dir)
    data_handler.mkdir(pep_dir)

    for i, (pep_id, coding_range) in enumerate(ref_peptides_dict.items(), 1):
        pep_id_safe = pep_id.replace("/", "_")
        nuc_path = nuc_dir / f"{i}_{pep_id_safe}"
        pep_path = pep_dir / f"{i}_{pep_id_safe}"
        data_handler.mkdir(nuc_path)
        data_handler.mkdir(pep_path)

        for sample, sample_seq in aligned_genome_seq_dict.items():
            fasta_file = f"{sample}_{pep_id_safe}.fa"
            nuc_seq = sample_seq[coding_range[0] - 1 : coding_range[1]]
            aa_seq = snp_analyzer.get_aa_seq(nuc_seq)

            with open(nuc_path / fasta_file, "w", encoding="utf-8") as n:
                n.write(f">{fasta_file}\n{nuc_seq}\n")

            with open(pep_path / fasta_file, "w", encoding="utf-8") as p:
                p.write(f">{fasta_file}\n{aa_seq}\n")

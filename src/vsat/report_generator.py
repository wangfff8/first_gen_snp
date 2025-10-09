#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import shutil
from pathlib import Path
import pandas as pd
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
    with open("snp_result.txt", "w", encoding="utf-8") as fp:
        for (sample, mutated_genome_seq), mutate_sites in zip(
            aligned_genome_seq_dict.items(), mutate_sites_list
        ):
            fp.write(f"{sample}:\n")

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
                        fp.write(templates[template_key].format(*snp_info))


            out_cds_sites = sorted(mutated_sites - in_cds_sites)
            for nuc_site in out_cds_sites:
                fp.write(templates[7].format(nuc_site, mutates_dict[nuc_site]))
            fp.write("\n")


def write_snp_results_to_xls(
    aligned_genome_seq_dict: dict[str, str],
    mutate_sites_list: list[str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_seq: str,
) -> None:
    """Writes a tab-separated SNP analysis report compatible with Excel."""
    with open("snp_result.xls", "w", encoding="utf-8") as fp:
        header = [
            "Sample",
            "Nuc Site",
            "Mutation",
            "Product",
            "AA Pos",
            "Ref Codon",
            "Mut Codon",
            "Ref AA",
            "Mut AA",
        ]
        fp.write("\t".join(header) + "\n")
        for (sample, mutated_genome_seq), mutate_sites in zip(
            aligned_genome_seq_dict.items(), mutate_sites_list
        ):
            if not mutate_sites:
                continue

            mutates_dict = {
                int(m.split("_")[0]): m.split("_")[1] for m in mutate_sites.split(",")
            }
            mutated_sites = set(mutates_dict.keys())
            in_cds_sites = set()

            for nuc_site in sorted(mutated_sites):
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
                        fp.write("\t".join(map(str, snp_info)) + "\n")
                        break

            out_cds_sites = sorted(mutated_sites - in_cds_sites)
            for nuc_site in out_cds_sites:
                fp.write(f"{sample}\t{nuc_site}\t{mutates_dict[nuc_site]}\n")
            fp.write("\n")

def align_stat(
    ref_peptides_dict: dict[str, tuple[int, int]],
    aligned_genome_seq_dict: dict[str, str],
    ref_name: str,
    ref_seq: str,
) -> None:
    """Calculates and writes alignment identity statistics to an Excel-compatible file."""
    with open("alignment_statistics.xls", "w", encoding="gbk") as stat:
        pep_features = [ref_name] + list(ref_peptides_dict.keys())
        stat.write("\t".join(h for h in pep_features) + "\n")

        for sample, sample_seq in aligned_genome_seq_dict.items():
            stat.write(f"{sample}\t{snp_analyzer.cal_identity(ref_seq, sample_seq)}\t")
            pep_align_rate = []
            for product_name, coding_range in ref_peptides_dict.items():
                ref_nuc_seq = ref_seq[coding_range[0] - 1 : coding_range[1]]
                qry_nuc_seq = sample_seq[coding_range[0] - 1 : coding_range[1]]
                ref_aa_seq = snp_analyzer.get_aa_seq(ref_nuc_seq)
                qry_aa_seq = snp_analyzer.get_aa_seq(qry_nuc_seq)
                pep_align_rate.append(snp_analyzer.cal_identity(ref_aa_seq, qry_aa_seq))
            stat.write("\t".join(pep_align_rate) + "\n")


def store_sequence(
    assembly_dir: str | Path,
    aligned_genome_seq_dict: dict[str, str],
    ref_peptides_dict: dict[str, tuple[int, int]],
    ref_genome_file: str | Path,
    out_seq_dir: str | Path,
) -> None:
    """
    Stores the processed sequences, categorized by product, into a directory structure.
    """
    assembly_dir = Path(assembly_dir)
    out_seq_dir = Path(out_seq_dir)
    ref_genome_file = Path(ref_genome_file)
    data_handler.mkdir(out_seq_dir)

    genome_dir = out_seq_dir / "genome"
    data_handler.mkdir(genome_dir)
    if ref_genome_file.exists():
        shutil.copy(ref_genome_file, genome_dir)

    for sample_seq_file in assembly_dir.glob("*.seq"):
        if sample_seq_file.exists():
            shutil.copy(sample_seq_file, genome_dir)

    nuc_dir = out_seq_dir / "gene"
    pep_dir = out_seq_dir / "protein"
    data_handler.mkdir(nuc_dir)
    data_handler.mkdir(pep_dir)

    for i, (product_name, coding_range) in enumerate(ref_peptides_dict.items(), 1):
        product_name_safe = product_name.replace("/", "_")
        nuc_path = nuc_dir / f"{i}_{product_name_safe}"
        pep_path = pep_dir / f"{i}_{product_name_safe}"
        data_handler.mkdir(nuc_path)
        data_handler.mkdir(pep_path)

        for sample, sample_seq in aligned_genome_seq_dict.items():
            seq_id = f"{sample}_{product_name_safe}"
            nuc_seq = sample_seq[coding_range[0] - 1 : coding_range[1]]
            aa_seq = snp_analyzer.get_aa_seq(nuc_seq)

            with open(nuc_path / f"{seq_id}_nuc.fa", "w", encoding="utf-8") as n:
                n.write(f">{sample} {product_name_safe}\n{nuc_seq}\n")

            with open(pep_path / f"{seq_id}_pep.fa", "w", encoding="utf-8") as p:
                p.write(f">{sample} {product_name_safe}\n{aa_seq}\n")


def generate_html_report(
    gene_mafft_file: Path | None,
    protein_mafft_file: Path | None,
    gene_dir_name: str | None,
    protein_dir_name: str | None,
) -> None:
    """
    Generates a single HTML report by reading the previously generated .xls files,
    converting them to HTML tables, and embedding links for pre-generated
    sequence alignments.
    """
    snp_result_xls = Path("snp_result.xls")
    align_stat_xls = Path("alignment_statistics.xls")

    # --- Alignment Statistics Table ---
    if not align_stat_xls.exists():
        print(f"Warning: {align_stat_xls} not found. Skipping its section in HTML report.")
        align_stat_html = f"<p>{align_stat_xls} not found.</p>"
    else:
        align_stat_df = pd.read_csv(
            align_stat_xls, sep="\t", encoding="gbk", engine="python"
        )
        align_stat_html = align_stat_df.to_html(
            index=False, classes="table table-bordered table-hover", justify="center"
        )

    # --- SNP Result Table ---
    if not snp_result_xls.exists():
        print(f"Warning: {snp_result_xls} not found. Skipping its section in HTML report.")
        snp_html = f"<p>{snp_result_xls} not found.</p>"
    else:
        snp_df = pd.read_csv(
            snp_result_xls, sep="\t", encoding="utf-8", engine="python"
        ).dropna(how="all")
        snp_html = snp_df.to_html(
            index=False,
            classes="table table-bordered table-hover",
            na_rep="",
            justify="center",
        )

    # --- Genome Alignment Section (Link only) ---
    genome_mafft_path = Path("processed_sequences/genome/genome_aligned_colorized.html")
    if genome_mafft_path.exists():
        genome_alignment_html = (
            f'<h3>Genome 比对</h3>'
            f'<a href="{genome_mafft_path}" class="file-link" target="_blank">查看完整的 Genome 比对报告</a>'
        )
    else:
        genome_alignment_html = "<h3>Genome 比对</h3><p>Genome 比对执行失败或未生成报告。</p>"


    # --- Gene Alignment Section (Embedded with iframe) ---
    gene_alignment_html = ""
    if gene_mafft_file and gene_mafft_file.exists() and gene_dir_name:
        gene_alignment_html = (
            f'<h3>Gene 比对 ({gene_dir_name})</h3>'
            f'<a href="{gene_mafft_file}" class="file-link" target="_blank">在新窗口中打开 {gene_dir_name} 比对报告</a>'
            f'<iframe src="{gene_mafft_file}" style="width: 100%; height: 400px; border: 1px solid #ddd; border-radius: 4px; background-color: #fff;"></iframe>'
        )
    else:
        gene_alignment_html = (
            '<h3>Gene 比对</h3>'
            '<p>未在报告中嵌入特定的 Gene 比对 (例如 "1_..."), 但所有比对文件都已生成。</p>'
            '<a href="processed_sequences/gene" class="file-link">浏览所有 Gene 比对文件</a>'
        )

    # --- Protein Alignment Section (Embedded with iframe) ---
    protein_alignment_html = ""
    if protein_mafft_file and protein_mafft_file.exists() and protein_dir_name:
        protein_alignment_html = (
            f'<h3>Protein 比对 ({protein_dir_name})</h3>'
            f'<a href="{protein_mafft_file}" class="file-link" target="_blank">在新窗口中打开 {protein_dir_name} 比对报告</a>'
            f'<iframe src="{protein_mafft_file}" style="width: 100%; height: 400px; border: 1px solid #ddd; border-radius: 4px; background-color: #fff;"></iframe>'
        )
    else:
        protein_alignment_html = (
            '<h3>Protein 比对</h3>'
            '<p>未在报告中嵌入特定的 Protein 比对 (例如 "1_..."), 但所有比对文件都已生成。</p>'
            '<a href="processed_sequences/protein" class="file-link">浏览所有 Protein 比对文件</a>'
        )

    # --- Final HTML Template ---
    html_template = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Viral SNP 分析报告</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; line-height: 1.6; color: #333; margin: 0; padding: 0; background-color: #f4f4f4; }}
        .container {{ max-width: 1200px; margin: 20px auto; padding: 20px; background-color: #fff; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }}
        h1, h2 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h1 {{ text-align: center; }}
        table {{ width: 100%; border-collapse: collapse; margin-bottom: 20px; font-size: 0.9em; }}
        th, td {{ padding: 12px 15px; border: 1px solid #ddd; text-align: left; }}
        thead th {{ background-color: #e9ecef; font-weight: bold; }}
        tbody tr:nth-of-type(even) {{ background-color: #f9f9f9; }}
        tbody tr:hover {{ background-color: #f1f1f1; }}
        a {{ color: #3498db; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .file-link {{ display: inline-block; margin-bottom: 15px; background-color: #3498db; color: white; padding: 8px 12px; border-radius: 4px; font-size: 0.9em; }}
        .file-link:hover {{ background-color: #2980b9; }}
        .table-container {{ overflow-x: auto; }}
        .alignment-container {{ font-family: "Courier New", Courier, monospace; white-space: pre; overflow-x: auto; background-color: #fdfdfd; padding: 15px; border: 1px solid #eee; border-radius: 4px; font-size: 0.85em; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Viral SNP 分析报告</h1>
        <h2>测序比对结果统计</h2>
        <a href="{align_stat_xls.name}" class="file-link">下载 alignment_statistics.xls</a>
        <div class="table-container">{align_stat_html}</div>
        <h2>SNP 结果</h2>
        <a href="{snp_result_xls.name}" class="file-link">下载 snp_result.xls</a>
        <div class="table-container">{snp_html}</div>
        <h2>多序列比对</h2>
        {genome_alignment_html}
        {gene_alignment_html}
        {protein_alignment_html}
    </div>
</body>
</html>
"""
    report_file = Path("viral_snp_report.html")
    with report_file.open("w", encoding="utf-8") as f:
        f.write(html_template)

    print(f"HTML report generated: {report_file}")

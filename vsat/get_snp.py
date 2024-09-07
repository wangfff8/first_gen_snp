#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
此脚本旨在检测一代测序结果组装后的序列与参考序列（病毒基因组）之间的单一位点变异（SNPs）。
它能够识别变异位点在基因组上的确切位置，确定变异位点所编码的氨基酸，指出变异位点在编码蛋白中的位置，描述碱基的突变信息，判断突变是否导致了氨基酸的改变。
此外，该脚本还能计算序列比对的一致性，并提取保存相关序列。
注意：此脚本专门设计用于处理单一位点突变，不支持插入（Insertions）或缺失（Deletions）等复杂变异。
'''
import os
import re
import sys
import json
import shutil
import platform
import argparse

# 检查操作系统
if platform.system() == "Windows":
    encoding = "gbk"  # Windows系统使用GBK编码
else:
    encoding = "utf-8"  # 其他系统使用UTF-8编码

# 密码子与氨基酸的对应关系字典
aa_dict = {
    'GCT': ['A', 'Ala', '丙氨酸'],
    'GCC': ['A', 'Ala', '丙氨酸'],
    'GCA': ['A', 'Ala', '丙氨酸'],
    'GCG': ['A', 'Ala', '丙氨酸'],
    'CGT': ['R', 'Arg', '精氨酸'],
    'CGC': ['R', 'Arg', '精氨酸'],
    'CGA': ['R', 'Arg', '精氨酸'],
    'CGG': ['R', 'Arg', '精氨酸'],
    'AGA': ['R', 'Arg', '精氨酸'],
    'AGG': ['R', 'Arg', '精氨酸'],
    'AAT': ['N', 'Asn', '天冬酰胺'],
    'AAC': ['N', 'Asn', '天冬酰胺'],
    'GAT': ['D', 'Asp', '天冬氨酸'],
    'GAC': ['D', 'Asp', '天冬氨酸'],
    'TGT': ['C', 'Cys', '半胱氨酸'],
    'TGC': ['C', 'Cys', '半胱氨酸'],
    'CAA': ['Q', 'Gln', '谷氨酰胺'],
    'CAG': ['Q', 'Gln', '谷氨酰胺'],
    'GAA': ['E', 'Glu', '谷氨酸'],
    'GAG': ['E', 'Glu', '谷氨酸'],
    'GGT': ['G', 'Gly', '甘氨酸'],
    'GGC': ['G', 'Gly', '甘氨酸'],
    'GGA': ['G', 'Gly', '甘氨酸'],
    'GGG': ['G', 'Gly', '甘氨酸'],
    'CAT': ['H', 'His', '组氨酸'],
    'CAC': ['H', 'His', '组氨酸'],
    'ATT': ['I', 'Ile', '异亮氨酸'],
    'ATC': ['I', 'Ile', '异亮氨酸'],
    'ATA': ['I', 'Ile', '异亮氨酸'],
    'TTA': ['L', 'Leu', '亮氨酸'],
    'TTG': ['L', 'Leu', '亮氨酸'],
    'CTT': ['L', 'Leu', '亮氨酸'],
    'CTC': ['L', 'Leu', '亮氨酸'],
    'CTA': ['L', 'Leu', '亮氨酸'],
    'CTG': ['L', 'Leu', '亮氨酸'],
    'AAA': ['K', 'Lys', '赖氨酸'],
    'AAG': ['K', 'Lys', '赖氨酸'],
    'ATG': ['M', 'Met', '甲硫氨酸'],
    'TTT': ['F', 'Phe', '苯丙氨酸'],
    'TTC': ['F', 'Phe', '苯丙氨酸'],
    'CCT': ['P', 'Pro', '脯氨酸'],
    'CCC': ['P', 'Pro', '脯氨酸'],
    'CCA': ['P', 'Pro', '脯氨酸'],
    'CCG': ['P', 'Pro', '脯氨酸'],
    'TCT': ['S', 'Ser', '丝氨酸'],
    'TCC': ['S', 'Ser', '丝氨酸'],
    'TCA': ['S', 'Ser', '丝氨酸'],
    'TCG': ['S', 'Ser', '丝氨酸'],
    'AGT': ['S', 'Ser', '丝氨酸'],
    'AGC': ['S', 'Ser', '丝氨酸'],
    'ACT': ['T', 'Thr', '苏氨酸'],
    'ACC': ['T', 'Thr', '苏氨酸'],
    'ACA': ['T', 'Thr', '苏氨酸'],
    'ACG': ['T', 'Thr', '苏氨酸'],
    'TGG': ['W', 'Trp', '色氨酸'],
    'TAT': ['Y', 'Tyr', '酪氨酸'],
    'TAC': ['Y', 'Tyr', '酪氨酸'],
    'GTT': ['V', 'Val', '缬氨酸'],
    'GTC': ['V', 'Val', '缬氨酸'],
    'GTA': ['V', 'Val', '缬氨酸'],
    'GTG': ['V', 'Val', '缬氨酸'],
    'TAA': ['*', '*', '终止密码子'],
    'TAG': ['*', '*', '终止密码子'],
    'TGA': ['*', '*', '终止密码子'],
}

# 报告模板
templates = {
    1: "第 {1} 个碱基 ({2}) 发生突变, 不影响 {3} 区第 {4} 个氨基酸的表达 ({7}); ",
    2: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 发生突变, 不影响 {3} 区第 {4} 个氨基酸的表达 ({7}); ",
    3: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 及第 {10} 碱基 ({11}) 发生突变, 不影响 {3} 区第 {4} 个氨基酸的表达 ({7}); ",
    4: "第 {1} 个碱基 ({2}) 发生突变, 造成 {3} 区第 {4} 个氨基酸由{7}变成{8}; ",
    5: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 发生突变, 造成 {3} 区第 {4} 个氨基酸由{7}变成{8}; ",
    6: "第 {1} 个碱基 ({2}) 和第 {9} 碱基 ({10}) 及第 {10} 碱基 ({11}) 发生突变, 造成 {3} 区第 {4} 个氨基酸由{7}变成{8}; ",
    7: "第 {} 个碱基 ({}) 发生突变; ",
}


# 参考基因组
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
GENOME_DIR = os.path.join(SCRIPT_DIR, '..', 'genome')
REF_GENOME_FILE = os.path.join(GENOME_DIR, 'ref_genome.json')


def load_reference_genomes_from_json(genome_file: str) -> dict:
    with open(genome_file, 'r', encoding='utf-8') as file:
        return json.load(file)


def print_reference_genomes(genome_dict: dict):
    for k, v in genome_dict.items():
        print(f"{k:<10}{','.join(v)}")
        

def parse_arguments() -> tuple[str, str, str]:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This script is used to study viral SNP using first generation sequencing."
    )
    parser.add_argument(
        '-g', '--locus', help='virus locus, see in the directory: {}'.format(GENOME_DIR))
    parser.add_argument(
        '-s', '--assembly', help="assembled directory(using DNAMAN), ending in .seq")
    parser.add_argument(
        '-l', '--list', action='store_true', help="list available reference genomes")

    args = parser.parse_args()

    if args.list:
        genome_dict = load_reference_genomes_from_json(REF_GENOME_FILE)
        print_reference_genomes(genome_dict)
        exit()

    # Now check if both -g and -s are provided
    if args.locus is None or args.assembly is None:
        parser.error("Please provide both -g/--locus and -s/--assembly options.")

    locus = args.locus.strip()
    assemble_dir = str(args.assembly).strip()
    genome_file = os.path.join(GENOME_DIR, f"{locus}.fasta")
    genome_ann_file = os.path.join(GENOME_DIR, f"{locus}.gff3")
    
    return assemble_dir, genome_file, genome_ann_file


def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def read_genome_sequence(seq_file: str) -> tuple:
    """
    读取基因组序列文件并返回序列ID和序列内容

    Args:
    seq_file (str): 包含基因组序列的文件路径

    Returns:
    tuple: 包含序列ID和转换为大写的基因组序列的元组
    """
    with open(seq_file, "r", encoding="utf-8") as f:
        seq_id = ""
        seq = ""
        for line in f:
            if line.startswith(">"):
                seq_id = line.strip().split()[0][1:]
            else:
                seq += line.strip().upper()
    return seq_id, seq


def calculate_position_and_pad_with_n(ref_seq: str, query_seq: str, n: int = 10) -> str:
    """
    计算样品序列在参考序列上的相对位置，并返回由N填充左右两段缺失序列的序列

    Args:
    ref_seq (str): 参考序列
    query_seq (str): 待比对的样品序列
    n (int): 初始匹配长度，默认为10

    Returns:
    str: 在待比对序列前后用"N"填充的序列
    """
    times = n
    index = -1
    for _ in range(times):
        start = query_seq[:n]
        index = ref_seq.find(start)
        if index != -1:
            break
    else:
        assert index != -1, f"错误：{query_seq[:n]}与参考基因组存在错配！！！"

    # 计算需要填充的"N"的数量
    padding_left = index
    padding_right = len(ref_seq) - (index + len(query_seq))

    return "N" * padding_left + query_seq + "N" * padding_right


def find_mutation_info(ref_seq: str, query_seq: str) -> str:
    """
    查找样品序列相较参考序列的突变位置和替换信息

    Args:
    ref_seq (str): 参考序列
    query_seq (str): 待比对的样品序列

    Returns:
    str: 返回包含突变位置和替换信息的字符串，格式为"position_ref>query,position_ref>query,..."
    """
    mutations = []
    for i in range(len(ref_seq)):
        if query_seq[i] == "N":
            continue
        if ref_seq[i] != query_seq[i]:
            mutation = f"{i + 1}_{ref_seq[i]}>{query_seq[i]}"
            mutations.append(mutation)
    return ",".join(mutations)


def cal_identity(seq1, seq2):
    cons_num = sum(1 for x, y in zip(seq1, seq2) if x == y)
    return "{:.2%}".format(cons_num / len(seq1))


def get_aa_seq(nuc_seq):
    if n := (len(nuc_seq) % 3):
        nuc_seq = nuc_seq[:-n]
    aa_seq = ""
    for j in range(0, len(nuc_seq), 3):
        aa_seq += aa_dict.get(nuc_seq[j:j + 3], ['*', '*', '*'])[0]
    return aa_seq


def get_genome_seq_dict(dir_):
    """ 读取拼接后的序列文件，存入一个字典"""
    genome_seq_dict = dict()  # {"ref_seq": "seq", "sample1": "seq1", ...}
    genome_seq_files = [os.path.join(dir_, file) for file in os.listdir(dir_) 
                        if file.endswith(".seq")]

    for genome_file in genome_seq_files:
        seq_name, seq = read_genome_sequence(genome_file)
        genome_seq_dict[seq_name] = seq

    return genome_seq_dict


def align_sequences_and_find_mutations(ref_seq: str, genome_seq_dict: dict) -> tuple:
    """
    对字典中的序列与参考序列进行比对，以'N'填充左右两端缺失的序列，并查找变异信息

    Args:
    genome_seq_dict (dict): 包含待比对序列的字典，键为序列名称，值为序列字符串

    Returns:
    tuple: 以'N'填充后的序列字典，突变位置和替换信息列表
    """

    aligned_genome_seq_dict = {}
    mutate_sites_list = []

    for sample_name, sample_sequence in genome_seq_dict.items():
        aligned_sequence = calculate_position_and_pad_with_n(ref_seq, sample_sequence)
        aligned_genome_seq_dict[sample_name] = aligned_sequence
        mutate_sites_list.append(
            find_mutation_info(ref_seq, aligned_sequence))

    return aligned_genome_seq_dict, mutate_sites_list



def parse_gff(gff3_file: str) -> dict:
    """
    解析GFF3文件并提取编码区域信息

    Args:
    gff3_file (str): GFF3文件的路径

    Returns:
    dict: 包含编码区域信息的字典，键为编码区域的产物名称，值为起始和终止位置的元组
    """
    coding_dict = {}
    with open(gff3_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            tmp_list = line.strip().split("\t")
            if len(tmp_list) != 9:
                continue

            # global STRAND
            # if tmp_list[2] == "region":
            #     STRAND = tmp_list[6]

            if match := re.search('product=(.*?)(?:;|$)', tmp_list[8]):
                coding_dict[match.group(1)] = (
                    int(tmp_list[3]), int(tmp_list[4]))

    print(coding_dict)
    return coding_dict



def get_snp_info(nuc_site, ref_seq, mutated_genome_seq, coding_range):
    """
    获取SNP信息，包括氨基酸位置和突变前后的氨基酸信息

    Args:
    nuc_site (int): 核苷酸在基因组上的位置
    ref_seq (str): 参考基因组序列
    mutated_genome_seq (str): 比对参考基因组并在两端用N填充缺失核苷酸之后的样品基因组序列
    coding_range (tuple): 某个编码序列在基因组上的编码区间范围

    Returns:
    tuple: 包含氨基酸位置、突变前后的氨基酸信息等内容的元组
    """

    cod_site = nuc_site - coding_range[0] + 1

    quo, rem = divmod(cod_site, 3)
    if rem != 0:
        aa_n = quo + 1  # aa_n: 第几个氨基酸, rem: 氨基酸的第几个核苷酸
    else:
        aa_n = quo
    if rem == 0:
        ori_aa = ref_seq[nuc_site - 3:nuc_site]
        after_aa = mutated_genome_seq[nuc_site - 3:nuc_site]
    elif rem == 1:
        ori_aa = ref_seq[nuc_site - 1:nuc_site + 2]
        after_aa = mutated_genome_seq[nuc_site - 1:nuc_site + 2]
    else:
        ori_aa = ref_seq[nuc_site - 2:nuc_site + 1]
        after_aa = mutated_genome_seq[nuc_site - 2:nuc_site + 1]

    ori_pep = aa_dict.get(ori_aa, ['*', '*', '*'])[2]  # 0: 单字母, 1, 三字母, 2 中文
    after_pep = aa_dict.get(after_aa, ['*', '*', '*'])[2]
    # (sample, nuc_site,  mutates_dict[nuc_site], product_name)
    snp_info = (aa_n, ori_aa, after_aa, ori_pep, after_pep)
    return snp_info


def write_snp_results_to_txt(aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq):
    # global ref_seq

    res = open("snp_result.txt", "w", encoding="utf-8")
    for aligned_genome_seq, mutate_sites in zip(aligned_genome_seq_dict.items(), mutate_sites_list):
        sample, mutated_genome_seq = aligned_genome_seq
        # print(f"{sample}:")
        res.write(f"{sample}:\n")

        mutates_dict = {}
        if mutate_sites == '': continue
        for mutate in mutate_sites.split(","):
            mutate = mutate.split("_")  # position_ref>query
            mutates_dict[int(mutate[0])] = mutate[1]

        in_cds = set()  # 编码区变异位点
        mutated_sites = mutates_dict.keys()

        for nuc_site in mutated_sites:
            for product_name, coding_range in ref_peptides_dict.items():
                sub_site = []
                if nuc_site in range(coding_range[0], coding_range[1] + 1):
                    in_cds.add(nuc_site)
                    sub_site.append(nuc_site)

#        for product_name, coding_range in ref_peptides_dict.items():
#            itv = range(coding_range[0], coding_range[1] + 1)
#            sub_itv_list = [list(itv)[i:i + 3] for i in range(0, len(itv), 3)]
#
#            for sub_itv in sub_itv_list:
#                sub_site = []
#                for nuc_site in mutated_sites:
#                    if nuc_site in sub_itv:
#                        in_cds.add(nuc_site)
#                        sub_site.append(nuc_site)
#
                #print(sub_site)
                if sub_site:

                    if len(sub_site) == 1:
                        snp_info = (sample, sub_site[0], mutates_dict[sub_site[0]], product_name) + \
                            get_snp_info(sub_site[0], ref_seq, mutated_genome_seq, coding_range)

                        if snp_info[-1] == snp_info[-2]:
                            # print(templates[1].format(*snp_info))
                            res.write(templates[1].format(*snp_info))
                        else:
                            # print(templates[4].format(*snp_info))
                            res.write(templates[4].format(*snp_info))

                    elif len(sub_site) == 2:
                        snp_info1 = (sample, sub_site[0], mutates_dict[sub_site[0]], product_name) + \
                            get_snp_info(sub_site[0], ref_seq, mutated_genome_seq, coding_range)
                        snp_info2 = (sample, sub_site[1], mutates_dict[sub_site[1]], product_name) + \
                            get_snp_info(sub_site[1], ref_seq, mutated_genome_seq, coding_range)

                        if snp_info1[-1] == snp_info1[-2]:
                            # print(templates[2].format(*(snp_info1 + snp_info2[1:3])))
                            res.write(templates[2].format(*(snp_info1 + snp_info2[1:3])))
                        else:
                            # print(templates[5].format(*(snp_info1 + snp_info2[1:3])))
                            res.write(templates[5].format(*(snp_info1 + snp_info2[1:3])))

                    elif len(sub_site) == 3:
                        snp_info1 = (sample, sub_site[0], mutates_dict[sub_site[0]], product_name) + \
                            get_snp_info(sub_site[0], ref_seq, mutated_genome_seq, coding_range)
                        snp_info2 = (sample, sub_site[1], mutates_dict[sub_site[1]], product_name) + \
                            get_snp_info(sub_site[1], ref_seq, mutated_genome_seq, coding_range)
                        snp_info3 = (sample, sub_site[2], mutates_dict[sub_site[2]], product_name) + \
                            get_snp_info(sub_site[2], ref_seq, mutated_genome_seq, coding_range)

                        if snp_info1[-1] == snp_info1[-2]:
                            # print(templates[3].format(*(snp_info1 + snp_info2[1:3] + snp_info3[1:3])))
                            res.write(templates[3].format(*(snp_info1 + snp_info2[1:3] + snp_info3[1:3])))
                        else:
                            # print(templates[6].format(*(snp_info1 + snp_info2[1:3] + snp_info3[1:3])))
                            res.write(templates[6].format(*(snp_info1 + snp_info2[1:3] + snp_info3[1:3])))

        out_cds = sorted(list(set(mutated_sites) - in_cds))  # 非编码区变异位点

        for nuc_site in out_cds:
            # print(templates[7].format(nuc_site, mutates_dict[nuc_site]))
            res.write(templates[7].format(nuc_site, mutates_dict[nuc_site]))
        res.write("\n")
    res.close()
    return


def write_snp_results_to_xls(aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq):

    content = open("snp_result.xls", "w", encoding=encoding)
    for aligned_genome_seq, mutate_sites in zip(
            aligned_genome_seq_dict.items(), mutate_sites_list):
        # print("mutate_sites:", mutate_sites)
        sample, mutated_genome_seq = aligned_genome_seq
        mutates_dict = {}

        if mutate_sites == "": continue
        for mutate in mutate_sites.split(","):
            mutate = mutate.split("_")
            mutates_dict[int(mutate[0])] = mutate[1]

        in_cds = set()  # 编码区变异位点
        mutated_sites = mutates_dict.keys()

        for nuc_site in mutated_sites:
            for product_name, coding_range in ref_peptides_dict.items():
                if nuc_site in range(coding_range[0], coding_range[1] + 1):
                    in_cds.add(nuc_site)
                    snp_info = (sample, nuc_site, mutates_dict[nuc_site], product_name) + \
                        get_snp_info(nuc_site, ref_seq, mutated_genome_seq, coding_range)
                    content.write("\t".join([str(i) for i in snp_info]) + '\n')

        # print(in_cds)
        out_cds = sorted(list(set(mutated_sites) - in_cds))  # 非编码区变异位点
        # print(out_cds)
        for nuc_site in out_cds:
            # print(sample, nuc_site, mutates_dict[nuc_site])
            content.write(f"{sample}\t{nuc_site}\t{mutates_dict[nuc_site]}\n")
        content.write("\n")
    content.close()
    return


def align_stat(ref_peptides_dict, aligned_genome_seq_dict, ref_name, ref_seq):
    # global ref_name, ref_seq
    stat = open("测序比对结果统计.xls", "w", encoding=encoding)

    pep_features = [ref_name] + list(ref_peptides_dict.keys())
    stat.write("\t".join(h for h in pep_features) + "\n")

    for sample, sample_seq in aligned_genome_seq_dict.items():
        stat.write(f"{sample}\t{cal_identity(ref_seq, sample_seq)}\t")
        pep_align_rate = []
        for pep_id, coding_range in ref_peptides_dict.items():
            ref_nuc_seq = ref_seq[coding_range[0] - 1:coding_range[1]]
            qry_nuc_seq = sample_seq[coding_range[0] - 1:coding_range[1]]
            ref_aa_seq = get_aa_seq(ref_nuc_seq)
            qry_aa_seq = get_aa_seq(qry_nuc_seq)
            pep_align_rate.append(cal_identity(ref_aa_seq, qry_aa_seq))
        stat.write("\t".join(pep_align_rate) + "\n")
    stat.close()
    return


def store_sequence(assemble_dir, aligned_genome_seq_dict, ref_peptides_dict, out_seq_dir = "序列拼接"):
    mkdir(out_seq_dir)

    genome_dir = os.path.join(out_seq_dir, "genome")
    mkdir(genome_dir)
    shutil.copy(ref_genome_file, genome_dir)
    sample_genome_files = [os.path.join(assemble_dir, file) for file in os.listdir(assemble_dir) if file.endswith(".seq")]
    for sample_seq_file in sample_genome_files:
        shutil.copy(sample_seq_file, genome_dir)

    nuc_dir = f"{out_seq_dir}/nucleotide"
    pep_dir = f"{out_seq_dir}/protein"
    mkdir(nuc_dir)
    mkdir(pep_dir)

    for sample, sample_seq in aligned_genome_seq_dict.items():
        i = 1

        for pep_id, coding_range in ref_peptides_dict.items():
            pep_id = str.replace(pep_id, "/", "_")
            nuc_path = f"{nuc_dir}/{str(i)}_{pep_id}"
            pep_path = f"{pep_dir}/{str(i)}_{pep_id}"
            mkdir(nuc_path)
            mkdir(pep_path)
            i += 1

            fasta_file = sample + "_" + pep_id
            with open(f"{nuc_path}/{fasta_file}.fa", "w") as n, \
                    open(f"{pep_path}/{fasta_file}.fa", "w") as p:

                nuc_seq = sample_seq[coding_range[0] - 1:coding_range[1]]
                n.write(">" + fasta_file + "\n" + nuc_seq + "\n")

                aa_seq = get_aa_seq(nuc_seq)
                p.write(">" + fasta_file + "\n" + aa_seq + "\n")

    return


def main():
    global ref_genome_file, aligned_genome_seq_dict, ref_peptides_dict

    assemble_dir, ref_genome_file, genome_ann_file = parse_arguments() 

    ref_peptides_dict = parse_gff(genome_ann_file)

    ref_name, ref_seq = read_genome_sequence(ref_genome_file)

    genome_seq_dict = get_genome_seq_dict(assemble_dir)

    aligned_genome_seq_dict, mutate_sites_list = align_sequences_and_find_mutations(ref_seq, genome_seq_dict)

    write_snp_results_to_txt(aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq)

    write_snp_results_to_xls(aligned_genome_seq_dict, mutate_sites_list, ref_peptides_dict, ref_seq)

    align_stat(ref_peptides_dict, aligned_genome_seq_dict, ref_name, ref_seq)

    aligned_genome_seq_dict[ref_name] = ref_seq

    store_sequence(assemble_dir, aligned_genome_seq_dict, ref_peptides_dict)

    return


if __name__ == "__main__":
    
    main()

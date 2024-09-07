#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import argparse

from utils import *


def parse_args():
    parser = argparse.ArgumentParser(description="")

    # Required Arguments
    parser.add_argument("--rawdata", required=True, help="URL to the genome FASTA file for ref-seq.", metavar='')
    parser.add_argument("--split_data", required=True, help="URL to the genome FASTA file for ref-seq.", metavar='')
    parser.add_argument("--assemble_dir", required=True, help="URL to the genome GTF file for ref-seq.", metavar='')
    parser.add_argument("--id_map", required=True, help="URL to the genome FASTA file for ref-seq.", metavar='')
    args = parser.parse_args()
    return args


def get_sample_map(id_map):
    id_map_dict = {}
    if os.path.isfile(id_map):
        with open(id_map, 'r') as f:
            for line in f:
                line_list = line.strip().split()
                id_map_dict[line_list[0]] = line_list[1]
        
    return id_map_dict


def run_cap3(assemble_dir, sample_list):
    with open("{}/run_cap3.sh".format(assemble_dir), "w") as o:
        for sample in sample_list:
            sample_dir = os.path.join(assemble_dir, sample)
            make_dir(sample_dir)
            o.write("cd {}\n".format(sample_dir))
            o.write("cap3 {}/{}.fasta\n".format(sample_dir, sample))
            o.write("seqkit split -i {}.fasta.cap.contigs\n\n".format(sample))
    os.system("bash {}/run_cap3.sh".format(assemble_dir))


if __name__ == "__main__":
    args = parse_args()
    id_map_dict = get_sample_map(args.id_map)
    sample_list = list(id_map_dict.values())
    print(sample_list)
    split_data(args.rawdata, args.split_data, id_map_dict)
    combine_seq2fasta(args.split_data, args.assemble_dir)
    run_cap3(args.assemble_dir, sample_list)

#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import re
import sys
import shutil
from itertools import chain


def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def split_data(path, out_path, id_map):
    seq_list = []
    ab1_list = []
    pdf_list = []

    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            file_fullname = os.path.join(dirpath, filename)
            if filename.endswith('.seq'):
                seq_list.append(file_fullname)
            elif filename.endswith('.ab1'):
                ab1_list.append(file_fullname)
            elif filename.endswith('.pdf'):
                pdf_list.append(file_fullname)

    if seq_list and ab1_list and len(seq_list) == len(ab1_list):
        for old, new in id_map.items():
            sample_out_path = os.path.join(out_path, new)
            sample_seq_out_path = os.path.join(sample_out_path, 'seq')
            sample_ab1_out_path = os.path.join(sample_out_path, 'ab1')
            make_dir(sample_seq_out_path)
            make_dir(sample_ab1_out_path)
            if pdf_list:
                sample_pdf_out_path = os.path.join(sample_out_path, 'pdf')
                make_dir(sample_pdf_out_path)

            for file_ in chain(seq_list, ab1_list, pdf_list):
                if re.search(old, os.path.basename(file_)):
                    if file_.endswith('.seq'):
                        shutil.copy(file_, sample_seq_out_path)
                    elif file_.endswith('.ab1'):
                        shutil.copy(file_, sample_ab1_out_path)
                    elif file_.endswith('.pdf'):
                        shutil.copy(file_, sample_pdf_out_path)

    else:
        print('Please check the input_dir path;'
                ' there are no *.seq and *.ab1 files under it,'
                ' or their numbers are different.')
        exit()


def convert_seq2fasta(path, out_path):
    for dirpath, dirnames, filenames in os.walk(path):
        if filenames and filenames[0].endswith('.seq'):
            sample_name = os.path.basename(os.path.split(dirpath)[0])
            sample_out_path = os.path.join(out_path, sample_name)
            make_dir(sample_out_path)
            for seq_file in filenames:
                seq_file_fullname = os.path.join(dirpath, seq_file)
                seq_name = seq_file.split('.')[0]
                out_seq_name = seq_name + ".fasta"
                out_seq_file = os.path.join(sample_out_path, out_seq_name)
                with open(seq_file_fullname) as r, open(out_seq_file, "w") as o:
                    o.write(">" + seq_name + "\n")
                    for line in r:
                        o.write(line)


def combine_seq2fasta(path, out_path):
    for dirpath, dirnames, filenames in os.walk(path):
        if filenames and filenames[0].endswith('.seq'):
            sample_name = os.path.basename(os.path.split(dirpath)[0])
            sample_out_path = os.path.join(out_path, sample_name)
            make_dir(sample_out_path)
            out_seq_name = sample_name + ".fasta"
            out_seq_file = os.path.join(sample_out_path, out_seq_name)
            with open(out_seq_file, "w") as o:
                for seq_file in filenames:
                    seq_file_fullname = os.path.join(dirpath, seq_file)
                    seq_name = seq_file.split('.')[0]
                    with open(seq_file_fullname) as r:
                        o.write(">" + seq_name + "\n")
                        for line in r:
                            o.write(line.strip() + "\n")



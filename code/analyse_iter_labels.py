#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to analyse iter_labels.csv
separate them in two MSA
Usage: python3 analyse_iter_labels.py inputfile_csv MSA_1.fasta MSA_2.fasta
inputfile: str, csv file

"""
# import statements
from sys import argv

import csv

# functions definitions
def parse_csv(filename):
    """
    function to parse csv file from ouroboros 'labels_iter.cav'

    :param filename: str
    :return:
    intera_seqs: set, it contains seqs index which are interacting sequences
    """
    intera_seqs = set()
    with open(filename, encoding='utf-8') as f:
        reader = csv.reader(f)
        index_row = -1
        for row in reader:
            index_row += 1
            seq_z = float(((row[0]).strip().split())[-1])
            if seq_z > 0.5:
                 intera_seqs.add(index_row)
    return intera_seqs

def parse_MSA_file(filename):
    """function to parse fasta file

    :param filename: str
    :return: names: list of str
    """
    lines = open(filename).readlines()
    names = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            line = line.split('/')[0]
            line = line.split('>')[1]
            names.append(line)
    return names

def count_Arath(filename):
    lines = open(filename).readlines()
    count = 0
    names = []
    index_row = -1
    Arath_row = []
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            index_row += 1
            line = line.split('/')[0]
            name = line.split(">")[1]
            line = line.split('_')[1]
            if line == "ARATH":
                count += 1
                Arath_row.append(index_row)
                names.append(name)
    return count, names, Arath_row

def main():
    intera_seq = parse_csv(argv[1])
    #MSA_a_names = parse_MSA_file(argv[2])
    # MSA_b_names = parse_MSA_file(argv[3])
    # intera_pairs_list=[]
    # # find all interaction seqs index
    # for index in intera_seq:
    #     # find interaction pairs name
    #     intera_pairs = {MSA_a_names[index], MSA_b_names[index]}
    #     intera_pairs_list.append(intera_pairs)
    #list of set
    # print(intera_pairs_list)
    ARATH_num, names, Arath_row = count_Arath(argv[2])
    print('Predicted interactions in total is: ' + str(len(intera_seq)))
    Arath_predict = set(intera_seq)&set(Arath_row)
    print('Predicted interactions in Arabidopsis is: ' + str(len(Arath_predict)))
    print('Total Arabidopsis is ' + str(ARATH_num))

if __name__ == '__main__':
    main()
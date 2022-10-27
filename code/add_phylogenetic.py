#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to add phylogenetic effect to sequences.
Usage: python3 add_phylogenetic.py filename outfile aminoacid_num group_num
filename : str, file containing sequences that we need
outfile: str, file containing sequences in fasta file
group_num:int number of group,
aminoacid_num: int, number of amino acid
"""
# import statements
from sys import argv
# functions definitions
global AA_TABLE
AA_TABLE = {1:'A', 2:'R', 3:'N', 4:'D', 5:'C', 6:'E', 7:'Q', 8:'G', 9:'H',
            10:'I', 11:'L', 12:'K', 13:'M',
           14:'F', 15:'P', 16:'S', 17:'T', 18:'W', 19:'Y', 20:'V'}

def parse_fasta_file(filename):
    """function to parse fasta file to only sequences list

    :param filename: str, name of fasta file
    :return: sequences: list, list of sequences
    """
    lines = open(filename).readlines()
    sequences = []
    new_seq = ''
    for line in lines:
        line=line.strip()
        if line.startswith('>'):
            sequences.append(new_seq)
            new_seq = ''
        else:
            new_seq += line
    sequences.append(new_seq)
    sequences = sequences[1:]
    return sequences

def add_phy(group_num, aa_num, sequences):
    """function to add phyogenetic effect to sequences

    :param group_num: int, number of subset groups
    :param aa_num: int, number of same amino acid in one group
    :param sequences: list, seuqences without phylogeny
    :return: new_seqs, list, sequences with phylogeny
    """
    group_num = int(group_num)
    aa_num = int(aa_num)
    protein_num = len(sequences)
    sub_num = int(protein_num/group_num)
    new_seqs=[]
    for index, seq in enumerate(sequences):
        for i in range(1, group_num):
            if sub_num*(i-1) <= index < sub_num * (i):
                seq = seq[:-aa_num] + AA_TABLE[i] * aa_num
                new_seqs.append(seq)
        if (sub_num * (group_num-1)) <= index:
            seq = seq[:-aa_num] + AA_TABLE[group_num] * aa_num
            new_seqs.append(seq)
    return new_seqs

def write_fasta_file(filename, all_seqs):
    f = open(filename, 'w')
    for index, seq in enumerate(all_seqs):
        f.write('>seq' + str(index + 1) + '\n')
        f.write(seq + '\n')
    f.close()

def main():
    sequences = parse_fasta_file(argv[1])
    last_int = int(len(sequences)* 0.5)
    intera_seqs = sequences[:last_int]
    non_seqs = sequences[last_int:]
    new_inter_seqs = add_phy(argv[4], argv[3], intera_seqs)
    new_non_seqs = add_phy(argv[4], argv[3], non_seqs)
    new_seqs = new_inter_seqs + new_non_seqs
    write_fasta_file(argv[2], new_seqs)

if __name__ == '__main__':
    main()
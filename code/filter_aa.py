#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to filter sequences containing amino acids with X.
Usage: python3 filter_aa.py filename outfile
filename : str, fasta file containing X AAs.
outfile: str, fasta file with no X AAs
"""

# import statements
from sys import argv
# functions definitions
def parse_fasta_file(filename):
    """function to parse fasta file to only sequences list

    :param filename: str, name of fasta file
    :return: sequences: list, list of sequences
    """
    lines = open(filename).readlines()
    names = []
    sequences = []
    new_seq = ''
    for line in lines:
        line=line.strip()
        if line.startswith('>'):
            names.append(line)
            sequences.append(new_seq)
            new_seq = ''
        else:
            new_seq += line
    sequences.append(new_seq)
    sequences = sequences[1:]
    return names, sequences

def delete_X_aa(names_list, sequence_list):
    """function to delete amino acids containing X

    :return: new_name_list: list of str,
    new_seq_list: list of str
    """
    new_name_list = []
    new_seq_list = []
    for index, seq in enumerate(sequence_list):
        if 'X' not in seq:
            new_seq_list.append(seq)
            new_name_list.append(names_list[index])
    return new_name_list, new_seq_list

def write_fasta_file(filename, seq_list, name_list):
    f = open(filename, 'w')
    for index, seq in enumerate(seq_list):
        f.write(name_list[index] + '\n')
        f.write(seq + '\n')
    f.close()

def main():
    names, sequences = parse_fasta_file(argv[1])
    new_name_list, new_seq_list = delete_X_aa(names, sequences)
    write_fasta_file(argv[2], new_seq_list, new_name_list)

if __name__ == '__main__':
    main()
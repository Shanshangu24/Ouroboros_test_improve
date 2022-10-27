#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to concatenate_alignments containing different domains
Usage: python3 concatenate_alignments.py filename
filename: str, fasta file from pfam
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
        line = line.strip()
        if line.startswith('>'):
            name = (line.split('/'))[0]
            names.append(name)
            sequences.append(new_seq)
            new_seq = ''
        else:
            new_seq += line
    sequences.append(new_seq)
    sequences = sequences[1:]
    return names, sequences

#def extra_seqs(name, names_list, seqs_list):
#indices = [i for i, x in enumerate(names_list) if x == name]


def main():
    names_1, seqs_1 = parse_fasta_file(argv[1])
    #print(names_1[:10])
    #print(seqs_1[:10])
    #print(len(set(names_1)))
    names_2, seqs_2 = parse_fasta_file(argv[2])
    #print(len(names_2))
    #print(len(set(names_2)))
    #names_3, seqs_3 = parse_fasta_file(argv[3])
    #print(len(names_3))
    #print(len(set(names_3)))
    same_name = set(names_1)&set(names_2)
    #print(same_name)
    f = open(argv[3], 'w')
    for name in same_name:
        new_seq = seqs_1[names_1.index(name)]+ seqs_2[names_2.index(name)]
        #print(names_1.index(name))
        f.write(name + '\n')
        f.write(new_seq + '\n')
    f.close()



if __name__ == '__main__':
    main()
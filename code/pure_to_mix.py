#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to make data from pure data filr to mix data file.
Usage python3 pure_to_mix.py filename1 filename2
filename1: str, a fasta file contains all interaction sequences
filename2: str, a fasta file contains non interaction sequences
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


def mix_sequences(all_sequence_num, int_fraction, int_sequences, non_int_sequences):
    """function to make certain num and certain fraction sequences

    :param all_sequence_num: int, num of all sequences
    :param int_fraction: float, fraction of interaction sequences
    :param int_sequences: list, list of interaction sequences
    :param non_int_sequences: list, list of non-interaction sequences
    :return: mix_int_seqs: list, list of interaction sequences that we need
    mix_nonint_seqs: list, list of non-interaction sequences that we need
    """
    int_sequence_num = int(all_sequence_num * int_fraction)
    nonint_sequence_num = all_sequence_num - int_sequence_num
    mix_int_seqs=int_sequences[:int_sequence_num]
    mix_nonint_seqs = non_int_sequences[:nonint_sequence_num]
    return mix_int_seqs, mix_nonint_seqs

def write_mix_file(filename, mix_int_seqs, mix_nonint_seqs):
    MSA = mix_int_seqs + mix_nonint_seqs
    f = open(filename, 'w')
    for index, seq in enumerate(MSA):
        f.write('>seq' + str(index + 1) + '\n')
        f.write(seq + '\n')
    f.close()

def main():
    """
    # 100, 4
    int_sequences_all_a = parse_fasta_file(argv[1])
    non_int_sequences_all_a = parse_fasta_file(argv[2])
    int_sequences_all_b = parse_fasta_file(argv[3])
    non_int_sequences_all_b = parse_fasta_file(argv[4])
    # 200, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(200, 0.5, int_sequences_all_a,
                  non_int_sequences_all_a)
    write_mix_file('mix_100_200_50%_4_a.fasta', mix_int_seqs_a, mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(200, 0.5,
                                                  int_sequences_all_b,
                                                  non_int_sequences_all_b)
    write_mix_file('mix_100_200_50%_4_b.fasta', mix_int_seqs_b, mix_nonint_seqs_b)
    # 500, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(500, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_500_50%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(500, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_500_50%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 1000, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(1000, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_1000_50%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(1000, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_1000_50%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 2000, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(2000, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_2000_50%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(2000, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_2000_50%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 3000, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(3000, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_3000_50%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(3000, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_3000_50%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 2000, 30%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(2000, 0.3,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_2000_30%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(2000, 0.3,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_2000_30%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 2000, 70%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(2000, 0.7,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_2000_70%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(2000, 0.7,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_2000_70%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 2000, 90%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(2000, 0.9,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_2000_90%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(2000, 0.9,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_2000_90%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)


    # 100, 8/15/30
    int_sequences_all_a = parse_fasta_file(argv[1])
    non_int_sequences_all_a = parse_fasta_file(argv[2])
    int_sequences_all_b = parse_fasta_file(argv[3])
    non_int_sequences_all_b = parse_fasta_file(argv[4])
    # 2000, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(2000, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_2000_50%_30_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(2000, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_2000_50%_30_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    """
    # 100, 4
    int_sequences_all_a = parse_fasta_file(argv[1])
    non_int_sequences_all_a = parse_fasta_file(argv[2])
    int_sequences_all_b = parse_fasta_file(argv[3])
    non_int_sequences_all_b = parse_fasta_file(argv[4])
    # 5000, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(5000, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_5000_50%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(5000, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_5000_50%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)
    # 10000, 50%
    mix_int_seqs_a, mix_nonint_seqs_a = mix_sequences(10000, 0.5,
                                                      int_sequences_all_a,
                                                      non_int_sequences_all_a)
    write_mix_file('mix_100_10000_50%_4_a.fasta', mix_int_seqs_a,
                   mix_nonint_seqs_a)
    mix_int_seqs_b, mix_nonint_seqs_b = mix_sequences(10000, 0.5,
                                                      int_sequences_all_b,
                                                      non_int_sequences_all_b)
    write_mix_file('mix_100_10000_50%_4_b.fasta', mix_int_seqs_b,
                   mix_nonint_seqs_b)


if __name__ == '__main__':
    main()
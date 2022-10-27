#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to creat synthetic data of Miguel literature.
Usage python3 create_pure_data.py
"""
# import statements

# functions definitions
import random
global AA
AA_LIST = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', \
           'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
AA_INT = ['AR', 'ND', 'CE', 'GH','IL', 'KM', 'FP', 'ST', 'WY', 'TV', 'PS', 'RT',
          'QV', 'EY', 'NS']


def create_interact_protein(interact_column, protein_length, protein_number):
    """function to create interaction proteins

    :param interact_column: list of int, it contains index of interaction
    residuals
    :param protein_length: int, length of protein
    :param protein_number: int, number of protein
    :return:
    MAS_1: list of str, they are proteins in MSA1
    MSA_2: list of str, they are proteins in MSA2
    interaction: dictionary, {1:['A','T']}
    """
    MSA_1 = []
    MSA_2 = []
    interaction = {}
    for col in interact_column:
        interact = random.choice(AA_INT)
        interaction[col] = [interact[0], interact[1]]
    for row in range(protein_number):
        seq_MSA_1 = ''
        seq_MSA_2 = ''
        for col in range(protein_length):
            if col in interact_column:
                seq = random.sample(interaction[col], 2)
                seq_MSA_1 += seq[0]
                seq_MSA_2 += seq[1]
            else:
                seq_MSA_1 += random.choice(AA_LIST)
                seq_MSA_2 += random.choice(AA_LIST)
        MSA_1.append(seq_MSA_1)
        MSA_2.append(seq_MSA_2)
    return MSA_1, MSA_2, interaction

def create_non_interact_protein(protein_length, protein_number, interaction, interact_column):
    """function to generate non-interaction proteins

    :param protein_length: int, length of protein
    :param protein_number: int, number of protein
    :param interaction:  dictionary, {1:['A','T']}
    :param interact_column: list of int, it contains index of interaction
    residuals
    :return:
    MAS_1: list of str, they are proteins in MSA1
    MSA_2: list of str, they are proteins in MSA2
    """
    MSA_1 = []
    MSA_2 = []
    for row in range(protein_number):
        seq_MSA_1 = ''
        seq_MSA_2 = ''
        int_pattern = 0
        for col in range(protein_length):
            if col in interact_column:
                seq_MSA_1 += random.choice(interaction[col])
                seq_MSA_2 += random.choice(interaction[col])
                if seq_MSA_1[col] != seq_MSA_2[col]:
                    int_pattern += 1
            else:
                seq_MSA_1 += random.choice(AA_LIST)
                seq_MSA_2 += random.choice(AA_LIST)
        if int_pattern <= len(interact_column)*0.5:
            MSA_1.append(seq_MSA_1)
            MSA_2.append(seq_MSA_2)
    return MSA_1, MSA_2

def write_in_fasta_file(MSA, filename):
    f = open(filename, 'w')
    for index, seq in enumerate(MSA):
        f.write('>seq' + str(index+1) + '\n')
        f.write(seq + '\n')
    f.close()


def main():
    """
    # 4 sites
    int_column = [2, 4, 6, 8]
    protein_length = 100
    protein_number = 10000
    MSA_1, MSA_2, interaction = create_interact_protein(int_column, protein_length, protein_number)
    protein_number = 20000
    n_MSA_1, n_MSA_2 = create_non_interact_protein(protein_length, protein_number, interaction, int_column)
    # print(MSA_1)
    # print(MSA_2)
    # print(n_MSA_1)
    # print(n_MSA_2)
    write_in_fasta_file(MSA_1, 'alnint_100_4_a.fasta')
    write_in_fasta_file(MSA_2, 'alnint_100_4_b.fasta')
    write_in_fasta_file(n_MSA_1, 'clean_nonint_100_4_a.fasta')
    write_in_fasta_file(n_MSA_2, 'clean_nonint_100_4_b.fasta')


    # 8 sites
    int_column = [2, 4, 6, 8, 10, 12, 14, 16]
    protein_length = 100
    protein_number = 10000
    MSA_1, MSA_2, interaction = create_interact_protein(int_column,
                                                        protein_length,
                                                        protein_number)
    protein_number = 20000
    n_MSA_1, n_MSA_2 = create_non_interact_protein(protein_length,
                                                   protein_number, interaction,
                                                   int_column)

    write_in_fasta_file(MSA_1, 'alnint_100_8_a.fasta')
    write_in_fasta_file(MSA_2, 'alnint_100_8_b.fasta')
    write_in_fasta_file(n_MSA_1, 'clean_nonint_100_8_a.fasta')
    write_in_fasta_file(n_MSA_2, 'clean_nonint_100_8_b.fasta')

    # 15 sites
    int_column = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    protein_length = 100
    protein_number = 10000
    MSA_1, MSA_2, interaction = create_interact_protein(int_column,
                                                        protein_length,
                                                        protein_number)
    protein_number = 20000
    n_MSA_1, n_MSA_2 = create_non_interact_protein(protein_length,
                                                   protein_number, interaction,
                                                   int_column)

    write_in_fasta_file(MSA_1, 'alnint_100_15_a.fasta')
    write_in_fasta_file(MSA_2, 'alnint_100_15_b.fasta')
    write_in_fasta_file(n_MSA_1, 'clean_nonint_100_15_a.fasta')
    write_in_fasta_file(n_MSA_2, 'clean_nonint_100_15_b.fasta')

    # 30 sites
    int_column = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32,
                  34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
    protein_length = 100
    protein_number = 10000
    MSA_1, MSA_2, interaction = create_interact_protein(int_column,
                                                        protein_length,
                                                        protein_number)
    protein_number = 20000
    n_MSA_1, n_MSA_2 = create_non_interact_protein(protein_length,
                                                   protein_number, interaction,
                                                   int_column)

    write_in_fasta_file(MSA_1, 'alnint_100_30_a.fasta')
    write_in_fasta_file(MSA_2, 'alnint_100_30_b.fasta')
    write_in_fasta_file(n_MSA_1, 'clean_nonint_100_30_a.fasta')
    write_in_fasta_file(n_MSA_2, 'clean_nonint_100_30_b.fasta')
    """

if __name__ == '__main__':
    main()
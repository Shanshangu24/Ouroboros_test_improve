#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to creat synthetic data of Miguel literature.
Usage python3 create_synthetic_data.py
"""
# import statements

# functions definitions
import random

global AA_TABLE
AA_TABLE = {'A': 0, 'R': 1, 'N': 2,
            'D': 3, 'C': 4,
            'E': 5, 'Q': 6, 'G': 7,
            'H': 8, 'I': 9, 'L': 10,
            'K': 11, 'M': 12, 'F': 13,
            'P': 14, 'S': 15, 'T': 16,
            'W': 17, 'Y': 18, 'V': 19,
            '-': 20}

global AA
AA_LIST = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', \
           'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
AA_INT = ['AR', 'ND', 'CE', 'GH','IL', 'KM', 'FP', 'ST', 'WY']


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
        for col in range(protein_length):
            if col in interact_column:
                stop = True
                while stop:
                    seq_1 = random.choice(AA_LIST)
                    seq_2 = random.choice(AA_LIST)
                    if (seq_1 not in interaction[col]) or (seq_2 not in interaction[col]):
                        seq_MSA_1 += seq_1
                        seq_MSA_2 += seq_2
                        stop = False
            else:
                seq_MSA_1 += random.choice(AA_LIST)
                seq_MSA_2 += random.choice(AA_LIST)
        MSA_1.append(seq_MSA_1)
        MSA_2.append(seq_MSA_2)
    return MSA_1, MSA_2


def write_in_fasta_file(MSA, MSA_N, filename):
    f = open(filename, 'w')
    all_seqs = MSA+MSA_N
    for index, seq in enumerate(all_seqs):
        f.write('>seq' + str(index+1) + '\n')
        f.write(seq + '\n')
    f.close()

def add_phylogenetic_effect(MSA_1, MSA_2):
    """function to add effect from 60-79.

    :param MSA_1:
    :param MSA_2:
    :return:
    """
    MSA_1_phy = []
    MSA_2_phy = []
    for seq in MSA_1:
        new_seq = seq[0:60]
        new_seq += 'A' * 20
        new_seq += seq[80:]
        MSA_1_phy.append(new_seq)
    for seq in MSA_2:
        new_seq = seq[0:60]
        new_seq += 'A' * 20
        new_seq += seq[80:]
        MSA_2_phy.append(new_seq)
    return MSA_1_phy, MSA_2_phy


def main():
    random.seed(1)
    """"
    # 1.total 100 len, total 200 num, 50%, (2,4,6,8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 200
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column, protein_length, protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length, n_protein_number, interaction, interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_200_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_200_50%_4_b.fasta')

    #1. total 100 len, total 500 num, 50%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 500
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_500_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_500_50%_4_b.fasta')

    # 1. total 100 len, total 1000 num, 50%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 1000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_1000_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_1000_50%_4_b.fasta')

    # 1. total 100 len, total 2000 num, 50%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 2000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_2000_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_2000_50%_4_b.fasta')

    # 1. total 100 len, total 3000 num, 50%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 3000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_3000_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_3000_50%_4_b.fasta')
    """
    """
    # 2. total 100 len, total 2000 num, 30%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 2000
    protein_number = int(total_protein_number * 0.3)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_2000_30%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_2000_30%_4_b.fasta')

    # 2. total 100 len, total 2000 num, 70%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 2000
    protein_number = int(total_protein_number * 0.7)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_2000_70%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_2000_70%_4_b.fasta')

    # 2. total 100 len, total 2000 num, 90%, (2, 4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 2000
    protein_number = int(total_protein_number * 0.9)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_2000_90%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_2000_90%_4_b.fasta')
    
    # 1.total 15 len, total 2000 num, 50%, (2,4)
    interact_column = [2, 4]
    protein_length = 15
    total_protein_number = 2000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_15_2000_50%_2_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_15_2000_50%_2_b.fasta')

# 1.total 15 len, total 4000 num, 50%, (2,4)
    interact_column = [2, 4]
    protein_length = 15
    total_protein_number = 4000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_15_4000_50%_2_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_15_4000_50%_2_b.fasta')
    """
    # 1.total 100 len, total 4000 num, 50%, (2,4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 5000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_5000_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_5000_50%_4_b.fasta')

    # 1.total 100 len, total 4000 num, 50%, (2,4, 6, 8)
    interact_column = [2, 4, 6, 8]
    protein_length = 100
    total_protein_number = 10000
    protein_number = int(total_protein_number * 0.5)
    n_protein_number = total_protein_number - protein_number
    MSA_1, MSA_2, interaction = create_interact_protein(interact_column,
                                                        protein_length,
                                                        protein_number)
    MSA_1_n, MSA_2_n = create_non_interact_protein(protein_length,
                                                   n_protein_number,
                                                   interaction,
                                                   interact_column)

    write_in_fasta_file(MSA_1, MSA_1_n, 'mix_100_10000_50%_4_a.fasta')
    write_in_fasta_file(MSA_2, MSA_2_n, 'mix_100_10000_50%_4_b.fasta')



if __name__ == "__main__":
    main()
#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to calculate similarities in MSAs
some functions from google
Usage: python3 calculate_similarities.py MSA_file threshold
inputfile:
MSA_file: MSA_exact
threshold: number, similarity score to remove sequences
"""

import numpy as np

from sys import argv

def parse_fasta_file(filename):
    """function to parse fasta file to only sequences list

    :param filename: str, name of fasta file
    :return: sequences: list, list of sequences
    names: list of names
    species: list of species name
    """
    lines = open(filename).readlines()
    names = []
    species = []
    sequences = []
    new_seq = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            try:
                name = line.split('>')[-1]
                specie_1 = line.split(' ')[0]
                specie_2 = specie_1.split('_')[1]
                specie = specie_2.split('/')[0]
                names.append(name)
                species.append(specie)
                sequences.append(new_seq)
                new_seq = ''
            except IndexError:
                continue
        else:
            new_seq += line
    sequences.append(new_seq)
    sequences = sequences[1:]
    return names, species, sequences


def calculate_similarity(seq_1, seq_2):
    numerator = 0
    denominator = 0
    for index , aa in enumerate(seq_1):
        aa_2 = seq_2[index]
        if (aa != '_') and (aa_2 != '_'):
            if aa == aa_2:
                numerator += 1
                denominator += 1
            else:
                denominator += 1
    similarity = numerator / denominator
    return similarity

# def sortbygroup(seqs_list, percent = 0.75):
#     groups = []
#     # go through every seq in sequences
#     for seq in seqs_list:
#         match = False
#         for g in range(len(groups)):
#             # find every group in groups
#             group = groups[g]
#             # find the first seq in groups
#             parent = group[0]
#             try:
#                 similarity = calculate_similarity(parent, seq)
#                 if similarity >= percent:
#                     group.append(seq)
#                     group.sort()
#                     match = True
#             except:
#                 pass
#         if not match:
#             groups.append([seq])
#     return groups


def buildSimilarityMatrix(seqs_lst, threshold ):
    numofSamples = len(seqs_lst)
    matrix = np.zeros(shape=(numofSamples, numofSamples))
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            similarity = calculate_similarity(seqs_lst[i], seqs_lst[j])
            if similarity > threshold:
                matrix[i, j] = 1
    return matrix

# function to calculate sum of row
def sumRow(matrix, i):
    return np.sum(matrix[i, :])

# function to find largest sum of row
def determineRow(matrix):
    maxNumOfOnes = -1
    row = -1
    for i in range(len(matrix)):
        if maxNumOfOnes < sumRow(matrix, i):
            maxNumOfOnes = sumRow(matrix, i)
            row = i
    return row


def addIntoGroup(matrix, ind):
    change = True
    indexes = []
    # add elements when it equals 1
    for col in range(len(matrix)):
        if matrix[ind, col] == 1:
            indexes.append(col)
    while change == True:
        change = False
        numIndexes = len(indexes)
        for i in indexes:
            for col in range(len(matrix)):
                if matrix[i, col] == 1:
                    if col not in indexes:
                        indexes.append(col)
        numIndexes2 = len(indexes)
        if numIndexes != numIndexes2:
            change = True
    return indexes


def deleteChosenRowsAndCols(matrix, indexes):
    for i in indexes:
        matrix[i, :] = 0
        matrix[:, i] = 0
    return matrix


def categorizeIntoClusters(matrix):
    groups = []
    while np.sum(matrix) > 0:
        row = determineRow(matrix)
        indexes = addIntoGroup(matrix, row)
        groups.append(indexes)
        matrix = deleteChosenRowsAndCols(matrix, indexes)
    return groups


def write_fasta_file(filename, seq_list, name_list):
    f = open(filename, 'w')
    for index, seq in enumerate(seq_list):
        #print(index)
        #print(name_list[index])
        f.write( '>' + name_list[index] + '\n')
        #print(seq)
        f.write(seq + '\n')
    f.close()


def main():
    # they use index to correspond each other
    names_list, species_list, sequences_list = parse_fasta_file("../real_data_pfam/proteasome_prot_A_N/proteasome_proteasome_A_N_exact.fasta")
    #threshold = float(argv[2])
    threshold = 0.95

    names_array = np.array(names_list)
    species_array = np.array(species_list)
    sequences_array = np.array(sequences_list)

    species = set(species_list)
    # remove similar elements
    new_names = []
    new_seqs = []
    for element in species:
        letter_list = np.where(species_array == element)
        # find index of species
        species_index = letter_list[0] #[2, 4, 6]
        if len(species_index) == 1:
            #print(element)
            names = names_array[species_index]
            new_names.extend(names)
            #print(names)
            seqs = sequences_array[species_index]
            new_seqs.extend(seqs)
            #print(seqs)
        else:
            if element == 'ARATH':
                names = names_array[species_index]
                print(names)
                seqs = sequences_array[species_index]
                #print(seqs)
                matrix = buildSimilarityMatrix(seqs, threshold)
                groups = categorizeIntoClusters(matrix)
                for group in groups:
                    # group is the index of groups
                    new_names.append(names[group[0]])
                    new_seqs.append(seqs[group[0]])
                print(element)
                print(groups)
    # print(len(new_names))
    # print(len(new_seqs))

    #write_fasta_file(argv[3], new_seqs, new_names)








if __name__ == '__main__':
    main()



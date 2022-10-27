#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to find sequences in the same species and
separate them in two MSA (within protein family)
Usage: python3 find_same_species.py filename MSA_1.fasta MSA_2.fasta
filename : str, fasta file.
MSA_1.fasta: str, fasta file with species corresponding to MSA_2
MSA_2.fasta: str, fasta file with species corresponding to MSA_1
"""

# import statements
from sys import argv
# functions definitions
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
        line=line.strip()
        if line.startswith('>'):
            names.append(line)
            specie_1 = line.split('/')[0]
            specie = specie_1.split('_')[1]
            species.append(specie)
            sequences.append(new_seq)
            new_seq = ''
        else:
            new_seq += line
    sequences.append(new_seq)
    sequences = sequences[1:]
    return names, species, sequences

def check_repeat_species(species_list):
    """function to sort sequence according to species

    :param species_list: list of species names
    :return:
    species_dict: dictionary contains species with multiple sequences
    species_dict[species_name] = [seq1_index, seq2_index, ...]
    """
    list1 = species_list
    # put the element in list2 if we find it at first time
    list2 = {}
    # put the element repeat many times in list3
    list3 = {}
    for index, element in enumerate(list1):
        if element not in list2:
            list2[element] = [index]
        # if the element repeat many times, we only store it once
        elif element in list2 and element not in list3:
            list3[element] = list2[element]
            list3[element].append(index)
        elif element in list2 and element in list3:
            list3[element].append(index)
    species_dict = list3
    return species_dict

def arrange_MSA_file_random(seqs_dict, seqs_list, name_list):
    """function to assign sequences into two MSA

    :param seqs_dict: dictionary contains species with multiple sequences
    species_dict[species_name] = [seq1_index, seq2_index, ...]
    :param seqs_list: list of sequences
    :param name_list: list of sequences names
    :return:
    MSA_1_name: list of MSA1 name
    MSA_1_seq: list of MSA1 sequences
    MSA_2_name: list of MSA2 name
    MSA_2_seq: list of MSA2 sequences
    """
    MSA_1_name = []
    MSA_1_seq = []
    MSA_2_name = []
    MSA_2_seq = []
    for species in seqs_dict:
        if len(seqs_dict[species]) % 2 == 0:
            for index_index, index in enumerate(seqs_dict[species]):
                if index_index % 2 == 0:
                    MSA_1_name.append(name_list[index])
                    MSA_1_seq.append(seqs_list[index])
                elif index_index % 2 != 0:
                    MSA_2_name.append(name_list[index])
                    MSA_2_seq.append(seqs_list[index])
        else:
            for index_index, index in enumerate(seqs_dict[species][:-1]):
                if index_index % 2 == 0:
                    MSA_1_name.append(name_list[index])
                    MSA_1_seq.append(seqs_list[index])
                elif index_index % 2 != 0:
                    MSA_2_name.append(name_list[index])
                    MSA_2_seq.append(seqs_list[index])
            for index_index, index in enumerate(seqs_dict[species][-2:]):
                if index_index % 2 == 0:
                    MSA_1_name.append(name_list[index])
                    MSA_1_seq.append(seqs_list[index])
                elif index_index % 2 != 0:
                    MSA_2_name.append(name_list[index])
                    MSA_2_seq.append(seqs_list[index])
    return MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq

def write_fasta_file(filename, seq_list, name_list):
    f = open(filename, 'w')
    for index, seq in enumerate(seq_list):
        f.write(name_list[index] + '\n')
        f.write(seq + '\n')
    f.close()


def main():
    names, species, sequences = parse_fasta_file(argv[1])
    species_dict = check_repeat_species(species)
    print(len(species_dict))
    print(len(species))
    print(species_dict)
    MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq = arrange_MSA_file_random(species_dict, sequences, names)
    write_fasta_file(argv[2], MSA_1_seq, MSA_1_name)
    write_fasta_file(argv[3], MSA_2_seq, MSA_2_name)

if __name__ == '__main__':
    main()
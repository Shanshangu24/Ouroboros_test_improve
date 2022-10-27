#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to find sequences in the same species and
separate them in two MSA (two protein family)
Usage: python3 find_same_species.py file1.fasta file2.fasta MSA_1.fasta MSA_2.fasta
file1 : str, fasta file.
file2: str. fasta file.
MSA_1.fasta: str, fasta file with species corresponding to MSA_2
MSA_2.fasta: str, fasta file with species corresponding to MSA_1

note: The file should correspond to the interactions columns
"""

# import statements
from sys import argv
import numpy as np
import random
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


def make_two_MSA(same_species, protein_1_species, protein_1_names, protein_1_seqs, protein_2_species, protein_2_names, protein_2_seqs):
    """function to make two file corresponding same species from two fasta file

    :param same_species: set. it contains the same species in two fasta file
    :param protein_1_species:list of str, it contains species name
    :param protein_1_names:list of str, it contains seqs' names
    :param protein_1_seqs:list of str, it contains seqs
    :param protein_2_species:the same as above three lines, but protein 2
    :param protein_2_names:
    :param protein_2_seqs:
    :return:
    MSA_1_name: list of MSA1 name
    MSA_1_seq: list of MSA1 sequences
    MSA_2_name: list of MSA2 name
    MSA_2_seq: list of MSA2 sequences

    note: this function contians all combination possibility.
    eg: for the species a, MSA_a has 3 seqs, MSA_b has 4 seqs. it will repeat like
    a1,   b1,
    a1,   b2,
    a1,   b3,
    a1,   b4
    ......
    """
    MSA_1_name = []
    MSA_1_seq = []
    MSA_2_name = []
    MSA_2_seq = []
    # convert list into numpy
    prot_2_name_array = np.array(protein_2_names)
    prot_2_seqs_array = np.array(protein_2_seqs)
    prot_1 = np.array(protein_1_species)
    prot_2 = np.array(protein_2_species)
    for species in same_species:
        prot_1_letter = np.where(prot_1 == species)
        prot_1_index = prot_1_letter[0] #[0,3,5]
        prot_2_letter = np.where(prot_2 == species)
        prot_2_index = prot_2_letter[0] #[0,2,4,6]
        for index in prot_1_index:
            MSA_1_name.extend([protein_1_names[index]]*len(prot_2_index))
            MSA_1_seq.extend([protein_1_seqs[index]]*len(prot_2_index))
            MSA_2_name.extend(prot_2_name_array[prot_2_index])
            MSA_2_seq.extend(prot_2_seqs_array[prot_2_index])
    return MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq

def make_two_MSA_random(same_species, protein_1_species, protein_1_names, protein_1_seqs, protein_2_species, protein_2_names, protein_2_seqs):
    """function to make two file corresponding same species from two fasta file

        :param same_species: set. it contains the same species in two fasta file
        :param protein_1_species:list of str, it contains species name
        :param protein_1_names:list of str, it contains seqs' names
        :param protein_1_seqs:list of str, it contains seqs
        :param protein_2_species:the same as above three lines, but protein 2
        :param protein_2_names:
        :param protein_2_seqs:
        :return:
        MSA_1_name: list of MSA1 name
        MSA_1_seq: list of MSA1 sequences
        MSA_2_name: list of MSA2 name
        MSA_2_seq: list of MSA2 sequences
        note: this function generates two MSA randomly
        for the species a, MSA_a has 3 seqs, MSA_b has 4 seqs. it will repeat like
        a1 b1
        a2 b2
        a3 b3
        a3 b1
        """
    MSA_1_name = []
    MSA_1_seq = []
    MSA_2_name = []
    MSA_2_seq = []
    # convert list into numpy
    prot_1 = np.array(protein_1_species)
    prot_2 = np.array(protein_2_species)
    for species in same_species:
        prot_1_letter = np.where(prot_1 == species)
        prot_1_index = prot_1_letter[0]  # [0,3,5]
        prot_2_letter = np.where(prot_2 == species)
        prot_2_index = prot_2_letter[0]  # [0,2,4,6]
        if len(prot_1_index) == len(prot_2_index):
            for i in prot_1_index:
                MSA_1_name.append(protein_1_names[i])
                MSA_1_seq.append(protein_1_seqs[i])
            for j in prot_2_index:
                MSA_2_name.append(protein_2_names[j])
                MSA_2_seq.append(protein_2_seqs[j])
        elif len(prot_1_index) < len(prot_2_index):
            for i, index in enumerate(prot_2_index):
                # use remainder to get the index of index for short list
                ii_1 = i % len(prot_1_index)
                index_1 = prot_1_index[ii_1]
                MSA_1_name.append(protein_1_names[index_1])
                MSA_1_seq.append(protein_1_seqs[index_1])
                MSA_2_name.append(protein_2_names[index])
                MSA_2_seq.append(protein_2_seqs[index])
        else:
            for i, index in enumerate(prot_1_index):
                # use remainder to get the index of index for short list
                ii_2 = i % len(prot_2_index)
                index_2 = prot_2_index[ii_2]
                MSA_1_name.append(protein_1_names[index])
                MSA_1_seq.append(protein_1_seqs[index])
                MSA_2_name.append(protein_2_names[index_2])
                MSA_2_seq.append(protein_2_seqs[index_2])
    return MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq

def add_real_interaction(filename):
    """function to change real interactions from file to list of lists

    :param filename:
    :return:
    """
    real_inters = []
    lines = open(filename).readlines()
    for line in lines:
        inter_a = line.split()[0]
        inter_b = line.split()[1]
        real_inters.append([inter_a, inter_b])
    return real_inters

def extract_real_interactions(real_inters, protein_1_species, protein_1_names, protein_1_seqs, protein_2_species, protein_2_names, protein_2_seqs):
    MSA_1_name = []
    MSA_1_seq = []
    MSA_2_name = []
    MSA_2_seq = []
    for inter in real_inters:
        inter_a, inter_b = inter
        if (inter_a in protein_1_names) and (inter_b in protein_2_names):
            index_a = protein_1_names.index(inter_a)
            index_b = protein_2_names.index(inter_b)
            MSA_1_name.append(protein_1_names.pop(index_a))
            MSA_1_seq.append(protein_1_seqs.pop(index_a))
            MSA_2_name.append(protein_2_names.pop(index_b))
            MSA_2_seq.append(protein_2_seqs.pop(index_b))
            protein_1_species.pop(index_a)
            protein_2_species.pop(index_b)
    return protein_1_species, protein_1_names, protein_1_seqs, protein_2_species, protein_2_names, protein_2_seqs, MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq


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
    # it aims to generate all combinations
    # names_1, species_1, sequences_1 = parse_fasta_file(argv[1])
    # names_2, species_2, sequences_2 = parse_fasta_file(argv[2])
    # same_species = set(species_1)&set(species_2)
    # MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq = make_two_MSA(same_species, species_1, names_1, sequences_1, species_2, names_2, sequences_2)
    # # print(len(MSA_1_name))
    # # print(len(MSA_2_name))
    # # print(len(MSA_1_seq))
    # # print(len(MSA_2_seq))
    # write_fasta_file(argv[3], MSA_1_seq, MSA_1_name)
    # write_fasta_file(argv[4], MSA_2_seq, MSA_2_name)

    # #it aims to generate combinations randomly
    # names_1, species_1, sequences_1 = parse_fasta_file(argv[1])
    # print('Protein A has ' + str(len(names_1)) + ' sequences.')
    # #print(species_1)
    # names_2, species_2, sequences_2 = parse_fasta_file(argv[2])
    # print('Protein B has ' + str(len(names_2)) + ' sequences.')
    # # file : real interactions
    # real_inters = add_real_interaction(argv[3])
    # species_1, names_1, sequences_1, species_2, names_2, sequences_2, real_MSA_1_name, real_MSA_1_seq, real_MSA_2_name, real_MSA_2_seq = extract_real_interactions(real_inters, species_1,  names_1, sequences_1, species_2, names_2, sequences_2)
    # print('Protein A has ' + str(len(species_1)) + ' unknown pairs.')
    # print(len(names_1))
    # print('Portein A has ' + str(len(real_MSA_1_name)) + ' real pairs.')
    # print('Portein B has ' + str(len(species_2)) + ' unknown pairs.')
    # print(len(names_2))
    # print('Portein B has ' + str(len(real_MSA_2_name)) + ' real pairs.')
    # same_species = set(species_1) & set(species_2)
    # print('Same species is ' + str(len(same_species)))
    # MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq = make_two_MSA_random(same_species,
    #                                                             species_1,
    #                                                             names_1,
    #                                                             sequences_1,
    #                                                             species_2,
    #                                                             names_2,
    #                                                             sequences_2)
    #
    # final_MSA_1_seq = real_MSA_1_seq + MSA_1_seq
    # final_MSA_1_name = real_MSA_1_name + MSA_1_name
    # final_MSA_2_seq = real_MSA_2_seq + MSA_2_seq
    # final_MSA_2_name = real_MSA_2_name + MSA_2_name
    #
    # print('MSA_a has ' + str(len(final_MSA_1_name)))
    # print('MSA_b has ' + str(len(final_MSA_2_name)))
    # print(len(final_MSA_1_seq))
    # print(len(final_MSA_2_seq))
    # # 10346 seqs
    # write_fasta_file(argv[4], final_MSA_1_seq, final_MSA_1_name)
    # write_fasta_file(argv[5], final_MSA_2_seq, final_MSA_2_name)

    # #it aims to generate all combinations
    # names_1, species_1, sequences_1 = parse_fasta_file(argv[1])
    # print('Protein A has ' + str(len(names_1)) + ' sequences.')
    # # print(species_1)
    # names_2, species_2, sequences_2 = parse_fasta_file(argv[2])
    # print('Protein B has ' + str(len(names_2)) + ' sequences.')
    # # file : real interactions
    # real_inters = add_real_interaction(argv[3])
    # species_1, names_1, sequences_1, species_2, names_2, sequences_2, real_MSA_1_name, real_MSA_1_seq, real_MSA_2_name, real_MSA_2_seq = extract_real_interactions(
    #     real_inters, species_1, names_1, sequences_1, species_2, names_2,
    #     sequences_2)
    # print('Protein A has ' + str(len(species_1)) + ' unknown pairs.')
    # print(len(names_1))
    # print('Portein A has ' + str(len(real_MSA_1_name)) + ' real pairs.')
    # print('Portein B has ' + str(len(species_2)) + ' unknown pairs.')
    # print(len(names_2))
    # print('Portein B has ' + str(len(real_MSA_2_name)) + ' real pairs.')
    # same_species = set(species_1) & set(species_2)
    # print('Same species is ' + str(len(same_species)))
    # MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq = make_two_MSA(
    #     same_species,
    #     species_1,
    #     names_1,
    #     sequences_1,
    #     species_2,
    #     names_2,
    #     sequences_2)
    #
    #
    #
    # final_MSA_1_seq = real_MSA_1_seq + MSA_1_seq
    # final_MSA_1_name = real_MSA_1_name + MSA_1_name
    # final_MSA_2_seq = real_MSA_2_seq + MSA_2_seq
    # final_MSA_2_name = real_MSA_2_name + MSA_2_name
    #
    # print('MSA_a has ' + str(len(final_MSA_1_name)))
    # print('MSA_b has ' + str(len(final_MSA_2_name)))
    # print(len(final_MSA_1_seq))
    # print(len(final_MSA_2_seq))
    # write_fasta_file(argv[4], final_MSA_1_seq, final_MSA_1_name)
    # write_fasta_file(argv[5], final_MSA_2_seq, final_MSA_2_name)

    # it aims to generate all combinations and then random choice to save runtime
    names_1, species_1, sequences_1 = parse_fasta_file(argv[1])
    print('Protein A has ' + str(len(names_1)) + ' sequences.')
    # print(species_1)
    names_2, species_2, sequences_2 = parse_fasta_file(argv[2])
    print('Protein B has ' + str(len(names_2)) + ' sequences.')
    # file : real interactions
    real_inters = add_real_interaction(argv[3])
    species_1, names_1, sequences_1, species_2, names_2, sequences_2, real_MSA_1_name, real_MSA_1_seq, real_MSA_2_name, real_MSA_2_seq = extract_real_interactions(
        real_inters, species_1, names_1, sequences_1, species_2, names_2,
        sequences_2)
    print('Protein A has ' + str(len(species_1)) + ' unknown pairs.')
    print(len(names_1))
    print('Portein A has ' + str(len(real_MSA_1_name)) + ' real pairs.')
    print('Portein B has ' + str(len(species_2)) + ' unknown pairs.')
    print(len(names_2))
    print('Portein B has ' + str(len(real_MSA_2_name)) + ' real pairs.')
    same_species = set(species_1) & set(species_2)
    print('Same species is ' + str(len(same_species)))
    MSA_1_name, MSA_1_seq, MSA_2_name, MSA_2_seq = make_two_MSA(
        same_species,
        species_1,
        names_1,
        sequences_1,
        species_2,
        names_2,
        sequences_2)


    # random find index and then decrease list length
    random.seed(0)
    idx = list(range(len(MSA_1_name)))
    idx_reduce = random.sample(idx, 20000)
    #print(idx_reduce)
    idx_reduce.sort()
    MSA_1_seq_reduce = []
    MSA_1_name_reduce = []
    MSA_2_seq_reduce = []
    MSA_2_name_reduce = []
    for i in idx_reduce:
        MSA_1_seq_reduce.append(MSA_1_seq[i])
        MSA_1_name_reduce.append(MSA_1_name[i])
        MSA_2_seq_reduce.append(MSA_2_seq[i])
        MSA_2_name_reduce.append(MSA_2_name[i])

    final_MSA_1_seq = real_MSA_1_seq + MSA_1_seq_reduce
    final_MSA_1_name = real_MSA_1_name + MSA_1_name_reduce
    final_MSA_2_seq = real_MSA_2_seq + MSA_2_seq_reduce
    final_MSA_2_name = real_MSA_2_name + MSA_2_name_reduce

    print('MSA_a has ' + str(len(final_MSA_1_name)))
    print('MSA_b has ' + str(len(final_MSA_2_name)))
    print(len(final_MSA_1_seq))
    print(len(final_MSA_2_seq))
    write_fasta_file(argv[4], final_MSA_1_seq, final_MSA_1_name)
    write_fasta_file(argv[5], final_MSA_2_seq, final_MSA_2_name)


if __name__ == '__main__':
    main()





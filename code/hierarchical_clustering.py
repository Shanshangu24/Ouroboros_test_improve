#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to calculate similarities in MSAs using hierarchical clustering
Usage: python3 hierarchical_clustering.py inputfile
inputfile:
MSA_file: MSA_exact
"""

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree, cut_tree
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

def buildSimilarityMatrix(seqs_lst):
    numofSamples = len(seqs_lst)
    matrix = np.zeros(shape=(numofSamples, numofSamples))
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            similarity = calculate_similarity(seqs_lst[i], seqs_lst[j])
            matrix[i, j] = similarity
    return matrix

def buildDistanceMatrix(matrix):
    numofSamples = len(matrix)
    distance_matrix = np.zeros(shape=(numofSamples, numofSamples))
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            distance = 1 - matrix[i, j]
            distance_matrix[i, j] = distance
    return distance_matrix

def buildHierarchicalCluster(distance_matrix):
    linkage_matrix = linkage(squareform(distance_matrix), method='single')
    return linkage_matrix

def cutHierarchicalCluster(linkage_matrix):
    cuttree = cut_tree(linkage_matrix, height=0.05)
    list_same = []
    for i in cuttree:
        seq_index = [x for x in range(len(cuttree)) if cuttree[x] == i]
        list_same.append([i[0], seq_index])
    # key is the index of clusters, values is the index of sequences
    dict_seq = dict(list_same)
    return cuttree, dict_seq


def to_newick(linkage_matrix: np.ndarray, labels) -> str:
    """
    Outputs linkage matrix as Newick style string
    Parameters
    ----------
    linkage_matrix : np.ndarray
        condensed distance matrix
    labels : list of str, optional
        leaf labels
    Returns
    -------
    newick : str
        linkage matrix in newick format tree
    Source:
    https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    tree = to_tree(linkage_matrix, rd=False)

    def get_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.2f%s" % (
                leaf_names[node.id],
                parentdist - node.dist,
                newick
            )
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = get_newick(
                node.get_left(),
                newick,
                node.dist,
                leaf_names
            )
            newick = get_newick(
                node.get_right(),
                ",%s" % (newick),
                node.dist,
                leaf_names
            )
            newick = "(%s" % (newick)
            return newick

    newick = get_newick(tree, "", tree.dist, labels)

    return newick

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
    names_list, species_list, sequences_list = parse_fasta_file('../real_data_pfam/proteasome_prot_A_N/proteasome_proteasome_A_N_exact.fasta')
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
        species_index = letter_list[0]  # [2, 4, 6]
        if len(species_index) == 1:
            #print(element)
            names = names_array[species_index]
            new_names.extend(names)
            #print(names)
            seqs = sequences_array[species_index]
            new_seqs.extend(seqs)
            #print(seqs)
        else:
            names = names_array[species_index]
            print(names)
            seqs = sequences_array[species_index]
            #print(seqs)
            similarity_matrix = buildSimilarityMatrix(seqs)
            distance_matrix = buildDistanceMatrix(similarity_matrix)
            linkage_matrix = buildHierarchicalCluster(distance_matrix)
            #print(linkage_matrix)
            labels = names
            newick = to_newick(linkage_matrix, labels)
            #print(newick)
            cuttree, dict_seq = cutHierarchicalCluster(linkage_matrix)
            print(cuttree)
            for cluster in dict_seq:
                seq_index = dict_seq[cluster][0]
                new_names.append(names[seq_index])
                new_seqs.append(seqs[seq_index])
    print( 'The number of names is ' + str(len(new_names)))
    print( 'The number of sequences is' + str(len(new_seqs)))

    #write_fasta_file(argv[2], new_seqs, new_names)






if __name__ == '__main__':
    main()







#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to calculate similarities in MSAs using hierarchical clustering
Usage: python3 reweight_sequences.py MSA.fasta
inputfile:
MSA_file: MSA.fasta file when it is a input file
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
    """
    lines = open(filename).readlines()
    names = []
    sequences = []
    new_seq = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            try:
                name = line.split('>')[-1]
                names.append(name)
                sequences.append(new_seq)
                new_seq = ''
            except IndexError:
                continue
        else:
            new_seq += line
    sequences.append(new_seq)
    sequences = sequences[1:]
    return names, sequences

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

def buildHierarchicalCluster(distance_matrix, md):
    linkage_matrix = linkage(squareform(distance_matrix), method=md)
    return linkage_matrix

def cutHierarchicalCluster(linkage_matrix, ht):
    cuttree = cut_tree(linkage_matrix, height=ht)
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

def assign_weight(cluster_dictionary):
    # generate a weight dictionary {index of seq: weight}
    weight_dictionary = {}
    for cluster in cluster_dictionary:
        number = len(cluster_dictionary[cluster])
        weight = 1/number
        for seq_index in cluster_dictionary[cluster]:
            weight_dictionary[seq_index] = weight
    return weight_dictionary




def main():
    # they use index to correspond each other
    names_list, sequences_list = parse_fasta_file("../simulated_data_miguel/mix_100_2000_50%_4_10_3_a.fasta")
    print(len(names_list))
    print(len(sequences_list))
    similarity_matrix = buildSimilarityMatrix(sequences_list)
    #print(similarity_matrix)
    distance_matrix = buildDistanceMatrix(similarity_matrix)
    print(distance_matrix)
    method = 'single'
    linkage_matrix = buildHierarchicalCluster(distance_matrix, method)
    labels = names_list
    newick = to_newick(linkage_matrix, labels)
    #print(newick)
    height = 0.78
    cuttree, dict_seq = cutHierarchicalCluster(linkage_matrix, height)
    #print(cuttree)
    print(dict_seq)
    print(len(dict_seq))
    weight_dictionary = assign_weight(dict_seq)
    print(weight_dictionary)


if __name__ == "__main__":
    main()
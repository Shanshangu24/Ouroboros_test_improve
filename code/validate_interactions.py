#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to validate predicts interactions and real interactions.
Usage: python3 validate_interactions.py

"""

# import statements
from sys import argv
import analyse_iter_labels
# functions definitions
def extract_physical_interactions(lines):
    """function to extract only physical interaction to get interator names

        :param lines: list of str, it contains every line in BioGrid
        :return: arabidopsis_lines: list of str
        """
    new_lines = []
    for line in lines:
        line = line.strip().split('\t')
        Exp_sys_type = line[12]
        if Exp_sys_type == 'physical':
            interactor_A = line[23]
            interactor_B = line[26]
            if interactor_A != '-' and interactor_B != '-':
                new_line = interactor_A + '\t' + interactor_B
                new_lines.append(new_line)
    return new_lines

def extract_arabidopsis(lines):
    """function to extract arabidopsis with only physical interaction

    :param lines: list of str, it contains every line in BioGrid
    :return: arabidopsis_lines: list of str, it contains only arabidopsis
    """
    new_lines = []
    for line in lines:
        line = line.strip().split('\t')
        ID_1 = line[35]
        ID_2 = line[36]
        Exp_sys_type = line[12]
        if ID_1 == 'Arabidopsis thaliana (Columbia)' and ID_2 == 'Arabidopsis thaliana (Columbia)' and Exp_sys_type == 'physical':
            interactor_A = line[23]
            interactor_B = line[26]
            if interactor_A != '-' and interactor_B != '-':
                new_line = interactor_A + '\t' + interactor_B
                new_lines.append(new_line)
    return new_lines

def readdom(file):
    """function to generate the pfam file into protein: domain_domain

        :param file: str, name of pfam file
        :return: newdict, dictionary {protein: domain_domain}
        """
    f = open(file, 'r')
    data = f.readlines()
    f.close()
    dict = {}
    for line in data[3:]:
        s = line.split('\t')[:-1]
        prot = s[0]
        dom = s[5]
        if prot in dict:
            tmp = dict[prot]
        else:
            tmp = []
        if tmp.count(dom) == 0:
            tmp.append(dom)
        tmp.sort()
        dict[prot] = tmp
    newdict = {}
    for protein in dict.keys():
        tmp = ""
        for domain in dict[protein]:
            tmp = tmp + domain + "_"
        newdict[protein] = tmp[:-1]
    return newdict


def readdomain(file):
    """function to generate the pfam file into protein: domain_domain

            :param file: str, name of pfam file
            :return: newdict, dictionary {protein: domain_domain}
    """
    f = open(file, 'r')
    data = f.readlines()
    protodomain={}
    for line in data[3:]:
        s = line.split('\t')
        prot = s[0]
        dom = s[5]
        name = s[-1]
        if prot in protodomain:
            if dom in protodomain[prot]:
                continue
            else:
                protodomain[prot].append(dom)
        else:
            protodomain[prot]=[name, dom]
    return protodomain





def map_protein_domains(interaction_list, prot_dict, domain_1_set, domain_2_set):
    interaction_pairs_list = []
    for interaction in interaction_list:
        interactor_a = interaction.split('\t')[0]
        interactor_b = interaction.split('\t')[1]
        if (interactor_a in prot_dict) and interactor_b in prot_dict:
            domains_a = set(prot_dict[interactor_a][1:])
            domains_b = set(prot_dict[interactor_b][1:])
            if (domains_a == domain_1_set) and (domains_b == domain_2_set):
                interaction_pairs = {prot_dict[interactor_a][0],
                                     prot_dict[interactor_b][0]}
                interaction_pairs_list.append(interaction_pairs)
                continue
            elif (domains_b == domain_1_set) and (domains_a == domain_2_set):
                interaction_pairs = {prot_dict[interactor_a][0],
                                     prot_dict[interactor_b][0]}
                interaction_pairs_list.append(interaction_pairs)
                continue
    return interaction_pairs_list






def map_int_domains(interaction_list, prot_dict, domain_1_set, domain_2_set):
    """function to find interactions containing certain domains.

    :param interaction_list:
    :param prot_dict:
    :param domain_1_set:
    :param domain_2_set:
    :return:
    """
    interaction_pairs_list = []
    for interactions in interaction_list:
        interactor_a_list = (interactions.split('\t')[0]).split('|')
        interactor_b_list = (interactions.split('\t')[1]).split('|')
        for trembl_a in interactor_a_list:
            for trembl_b in interactor_b_list:
                if (trembl_a in prot_dict) and (trembl_b in prot_dict):
                    domains_a = set(prot_dict[trembl_a])
                    domains_b = set(prot_dict[trembl_b])
                    # domain 1 in domains a & domain 2 in domains b
                    same_domain_a = domain_1_set & domains_a
                    same_domain_b = domain_2_set & domains_b
                    if (len(same_domain_a) == len(domain_1_set)) and (len(same_domain_b) == len(domain_2_set)):
                        interaction_pairs = {prot_dict[trembl_a][0], prot_dict[trembl_b][0]}
                        interaction_pairs_list.append(interaction_pairs)
                        continue
                    # domain 1 in domains b & domain 2 in domains a
                    same_domain_1 = domain_1_set & domains_b
                    same_domain_2 = domain_2_set & domains_a
                    if (len(same_domain_1) == len(domain_1_set)) and (len(same_domain_2) == len(domain_2_set)):
                        interaction_pairs = {prot_dict[trembl_b][0], prot_dict[trembl_a][0]}
                        interaction_pairs_list.append(interaction_pairs)
                        continue
    return interaction_pairs_list



def main():

    lines = open("../real_interactions/BIOGRID-ALL-4.4.211.tab3.txt").readlines()
    lines = lines[1:]
    #physical_intera = extract_physical_interactions(lines)
    #print(physical_intera[:50])
    ara_physical_intera = extract_arabidopsis(lines)
    #print(ara_physical_intera)
    domain_dict = readdomain("../real_interactions/Arath_domain_map.tsv")
    #print(domain_dict)


    # aux (PF02309) and aux (PF02309) b3 (PF02362) auxin_resp (PF06507)
    domain_1_set = {"PF02309"}
    #domain_2_set = {"PF02309"}
    domain_2_set = {"PF02309", "PF02362", "PF06507"}
    interaction_pairs_list = map_int_domains(ara_physical_intera, domain_dict, domain_1_set, domain_2_set)
    #print(interaction_pairs_list)


    """
    # f-box (PF00646)  and skp1(PF01466) skp1_poz (PF03931)
    domain_1_set = {"PF00646"}
    domain_2_set = {"PF01466", "PF03931"}
    interaction_pairs_list = map_int_domains(ara_physical_intera, domain_dict,
                                             domain_1_set, domain_2_set)
    #print(interaction_pairs_list)
    """
    """
    # f-box (PF00646) FBA_1 (PF07734) and skp1(PF01466) skp1_poz (PF03931)
    domain_1_set = {"PF00646", "PF07734"}
    domain_2_set = {"PF01466", "PF03931"}
    interaction_pairs_list = map_int_domains(ara_physical_intera, domain_dict,
                                             domain_1_set, domain_2_set)
    """

    # predict interactions
    intera_seq = analyse_iter_labels.parse_csv(argv[1])
    MSA_a_names = analyse_iter_labels.parse_MSA_file(argv[2])
    MSA_b_names = analyse_iter_labels.parse_MSA_file(argv[3])
    ARATH_num, ARA_names = analyse_iter_labels.count_Arath(argv[2])
    intera_pairs_list = []
    # find all interaction seqs index
    for index in intera_seq:
        # find interaction pairs name
        intera_pairs = {MSA_a_names[index], MSA_b_names[index]}
        intera_pairs_list.append(intera_pairs)
    # list of set
    #print(intera_pairs_list)
    #print(ARA_names)

    # calculate percentage
    count = 0
    for intr in interaction_pairs_list:
        if intr in intera_pairs_list:
            count += 1
    #how many real and predicted interactions in ARATH
    print(count)
    #how many real interactions in ARATH
    print(len(interaction_pairs_list))
    #how many predicted interactions
    print(len(intera_pairs_list))
    #how many predicted interactions in ARATH
    print(ARATH_num)
    ara_per = count / len(interaction_pairs_list)
    pred_per = count / len(intera_pairs_list)
    print(ara_per)
    print(pred_per)


    """
    # f-box (PF00646)  and skp1(PF01466) skp1_poz (PF03931)
    domain_1_set = {"PF00646"}
    domain_2_set = {"PF01466", "PF03931"}
    interaction_pairs_list = map_int_domains(ara_physical_intera, domain_dict,
                                             domain_1_set, domain_2_set)
    print(interaction_pairs_list)

    # predict interactions
    intera_seq = analyse_iter_labels.parse_csv(argv[1])
    MSA_a_names = analyse_iter_labels.parse_MSA_file(argv[2])
    MSA_b_names = analyse_iter_labels.parse_MSA_file(argv[3])
    intera_pairs_list = []
    # find all interaction seqs index
    for index in intera_seq:
        # find interaction pairs name
        intera_pairs = {MSA_a_names[index], MSA_b_names[index]}
        intera_pairs_list.append(intera_pairs)

    # calculate percentage
    count = 0
    for intr in interaction_pairs_list:
        if intr in intera_pairs_list:
            count += 1
    ara_per = count / len(interaction_pairs_list)
    pred_per = count / len(intera_pairs_list)
    print(ara_per)
    print(pred_per)
    """

if __name__ == '__main__':
    main()
#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to find real interactions containing certain domains
Usage: python3 find_real_interactions.py

"""

# import statements
from sys import argv
# functions definitions
def extract_arabidopsis(filename):
    """function to extract arabidopsis with only physical interaction

    :param lines: list of str, it contains every line in BioGrid
    :return: arabidopsis_lines: list of str, it contains only arabidopsis
    """
    lines = open(filename).readlines()
    lines = lines[1:]
    new_lines = []
    for line in lines:
        line = line.strip().split('\t')
        ID_1 = line[35]
        ID_2 = line[36]
        Exp_sys_type = line[12]
        if (ID_1 == 'Arabidopsis thaliana (Columbia)') and (ID_2 == 'Arabidopsis thaliana (Columbia)') and (Exp_sys_type == 'physical'):
            interactor_A = line[23]
            interactor_B = line[26]
            if interactor_A != '-' and interactor_B != '-':
                new_line = (interactor_A, interactor_B)
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
    prot2name = {}
    for line in data[3:]:
        s = line.strip().split('\t')
        prot = s[0]
        dom = s[5]
        name = s[-1]
        prot2name[prot] = name
        if prot in dict:
            tmp = dict[prot]
        else:
            tmp = []
        if tmp.count(dom) == 0:
            tmp.append(dom)
        tmp.sort()
        dict[prot] = tmp
    newdict = {}
    for prot in dict.keys():
        tmp = ""
        for domain in dict[prot]:
            tmp = tmp + domain + "_"
        newdict[prot] = tmp[:-1]
    return newdict, prot2name

def select_certain_domain(prot2domain, prot2name, inter_list, domain1, domain2):
    domain1.sort()
    domain2.sort()
    tmp = ''
    pairs = []
    for domain in domain1:
        tmp = tmp + domain + '_'
    domain1 = tmp[:-1]
    tmp = ""
    for domain in domain2:
        tmp = tmp + domain + '_'
    domain2 = tmp[:-1]
    if domain1 < domain2:
        target_domain = (domain1, domain2)
    else:
        target_domain = (domain2, domain1)
    print(target_domain)
    for i in inter_list:
        a, b = i
        if (a in prot2domain) and (b in prot2domain):
            dom1 = prot2domain[a]
            dom2 = prot2domain[b]
            if dom1 < dom2:
                tmp = (dom1, dom2)
                if tmp == target_domain:
                    # name for dom1, name for dom 2
                    inter_pairs = (prot2name[a], prot2name[b])
                    pairs.append(inter_pairs)
            else:
                tmp = (dom2, dom1)
                if tmp == target_domain:

                    inter_pairs = (prot2name[b], prot2name[a])
                    pairs.append(inter_pairs)
    pairs = set(pairs)
    return pairs

def write_pairs_file(pairs_set, filename ):
    with open(filename, 'w') as f:
        for pairs in pairs_set:
            a, b = pairs
            f.write(a + '\t' + b + '\n' )

def main():
    # read "Arath_domain_map.tsv"
    prot2dom, prot2name = readdom(argv[1])
    inter = extract_arabidopsis(argv[2])
    # print(len(inter))
    # print(inter[:10])

    # # aux (PF02309) and aux (PF02309) b3 (PF02362) auxin_resp (PF06507)
    # domain1 = ['PF02309']
    # domain2 = ['PF02309', 'PF02362', 'PF06507']
    # pairs = select_certain_domain(prot2dom, prot2name, inter, domain1, domain2)
    # print(len(pairs))
    # write_pairs_file(pairs, 'real_pairs_aux_aux_b3_auxin.txt')

    # # f-box (PF00646)  and skp1(PF01466) skp1_poz (PF03931)
    # domain1 = ['PF00646']
    # domain2 = ['PF01466', 'PF03931']
    # pairs = select_certain_domain(prot2dom, prot2name, inter, domain1, domain2)
    # print(len(pairs))
    # write_pairs_file(pairs, 'real_pairs_f_box_skp1_skp1_poz.txt')

    # # f-box (PF00646) FBA_1 (PF07734) and skp1(PF01466) skp1_poz (PF03931)
    # domain1 = ['PF00646', 'PF07734']
    # domain2 = ['PF01466', 'PF03931']
    # pairs = select_certain_domain(prot2dom, prot2name, inter, domain1, domain2)
    # print(len(pairs))
    # write_pairs_file(pairs, 'real_pairs_f_box_fba_1_skp1_skp1_poz.txt')

    # Proteasome(PF00227) and Proteasome(PF00227) Proteasome_A_N(PF10584)
    domain1 = ['PF00227']
    domain2 = ['PF00227', 'PF10584']
    pairs = select_certain_domain(prot2dom, prot2name, inter, domain1, domain2)
    print(len(pairs))
    write_pairs_file(pairs, 'real_pairs_Proteasome_Proteasome_A_N.txt')



if __name__ == '__main__':
    main()
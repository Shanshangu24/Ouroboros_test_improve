#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to find sequences with certain architecture
 map smart file to pfam file
Usage: python3 select_architecture.py

"""

# import statements
from sys import argv
# functions definitions
def extract_smart(filename):
    """function to extract names from smart

    :param filename:
    :return:
    notes: it means it only contains certain domains
    """
    real_names = set()
    lines = open(filename).readlines()
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            entry_name = ((line.split('/')[0]).split()[0]).split('|')[-1]
            real_names.add(entry_name)
    return real_names

def extract_domain_pfam(filename, real_names, new_file):
    lines = open(filename).readlines()
    fasta_dict = {}
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            name = (line.split('/')[0]).split('>')[-1]
            #print(name)
            fasta_dict[name] = ''
        else:
            fasta_dict[name] += line.replace('\n', '')
    count = 0
    with open(new_file, 'w') as f:
        for real_name in real_names:
            if real_name in fasta_dict:
                count += 1
                f.write('>' + real_name + '\n')
                f.write(fasta_dict[real_name] + '\n')
    return count


def main():
    real_names = extract_smart(argv[1])
    #print(real_names)
    count_a = extract_domain_pfam(argv[2], real_names, argv[3])
    print(count_a)

if __name__ == '__main__':
    main()


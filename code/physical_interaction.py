#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to find physical interactions of arabidopsis
thaliana from BioGrid tab file.
usage:
python3 physical_interaction.py inputfile
inputfile: str, name of input file from BioGrid
"""
# import statements
from sys import argv

# functions definitions
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

def write_file(lines, filename):
    lines.sort()
    with open(filename, 'w') as f:
        for line in lines:
            line = line.strip().split('\t')
            IA = line[0]
            IB = line[1]
            f.write(IA + '\t' + IB + '\n')

def extract_physical_interactions(lines):
    """function to extract only physical interaction to get doamin names

        :param lines: list of str, it contains every line in BioGrid
        :return: arabidopsis_lines: list of str
        """
    new_lines = []
    for line in lines:
        line = line.strip().split('\t')
        Exp_sys_type = line[12]
        if Exp_sys_type == 'physical':
            interactor_A = line[9]
            interactor_B = line[10]
            if interactor_A != '-' and interactor_B != '-':
                new_line = interactor_A + '\t' + interactor_B
                new_lines.append(new_line)
    return new_lines

def main():
    lines = open(argv[1]).readlines()
    lines = lines[1:]
    # new_lines = extract_arabidopsis(lines)
    # write_file(new_lines, 'Arabidopsis_int.txt')


if __name__ == '__main__':
    main()
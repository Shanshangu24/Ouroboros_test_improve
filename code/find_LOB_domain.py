#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to find physical interactions of LOB domains
usage:
python3 find_LOB_domain.py inputfile
inputfile: str, name of input file containing two columns of domain names
"""
# import statements
from sys import argv

# functions definitions
def parse_and_find_LOB(filename):
    lob_lob = 0
    lob_ints = []
    lines = open(filename).readlines()
    for line in lines:
        line = line.split('\t')
        domain_1 = line[0]
        domain_2 = line[1]
        if ('LATERAL ORGAN BOUNDARIES' in domain_1) or ('LOB' in domain_1) or ('lateral organ boundaries' in domain_1):
            count_domain_1 = 'LOB'
        else:
            count_domain_1 = domain_1
        if ('LATERAL ORGAN BOUNDARIES' in domain_2) or ('LOB' in domain_2) or ('lateral organ boundaries' in domain_2):
            count_domain_2 = 'LOB'
        else:
            count_domain_2 = domain_2
        if (count_domain_2 == 'LOB') or (count_domain_1 == 'LOB'):
            lob_int = [count_domain_1, count_domain_2]
            lob_ints.append(lob_int)
        if (count_domain_1 == 'LOB') and (count_domain_2 == 'LOB'):
            lob_lob += 1
    return lob_lob, lob_ints

def main():
    lob_lob, lob_ints = parse_and_find_LOB(argv[1])
    print(lob_lob)
    print(len(lob_ints))
    #print(lob_ints[:100])

if __name__ == '__main__':
    main()





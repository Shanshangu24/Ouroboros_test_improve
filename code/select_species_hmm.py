#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to extract species name from hmmsearch output file

Usage: python3 select_species_hmm.py
hmm_file1: str,  hmmsearch output file
hmm_file2: str, hmmsearch output file
"""

# import statements
from sys import argv
# functions definitions
def parse_hmmsearch(filename):
    lines = open(filename).readlines()
    timer = 0
    for line in lines:
        timer += 1
        if timer >=17:
            line = line.strip().split()
            seq_name =  line[8]
            


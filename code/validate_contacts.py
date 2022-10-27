#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to validate predicts interactions and real interactions.
Usage: python3 validate_interactions.py

"""

# import statements
from sys import argv
import pandas as pd
# functions definations
def parse_csv(filename):
    contact_matrix_df = pd.read_csv(filename)
    contact_matrix_df[contact_matrix_df != 0] = 1
    return contact_matrix_df






def main():
    contact_matrix_df = parse_csv(argv[1])
    print(contact_matrix_df)
    real_contact_matrix = pd.read_csv(argv[2])
    #contact_matrix_df.compare(real_contact_matrix)



if __name__ == '__main__':
    main()
#!/usr/bin/env python3

"""
Author；Shanshan GU
Description: This is a script to convert entry ID into entry names
Usage: python3 convert_entryname.py domainfile
"""
# import statements
import csv
from sys import argv
import pandas as pd
# functions definitions
def parse_pfam_domain(domainfile, outputfile):
    """function to read entry ID

        :param file:
        :return:
        """
    lines = open(domainfile).readlines()
    f = open(outputfile, 'w')
    for line in lines[3:]:
        entry_id = line.strip().split()[0]
        f.write(entry_id + "\t")
    f.close()




def main():
    #parse_pfam_domain(argv[1], "arath_entry_id.txt")
    domain_df = pd.read_table('../real_interactions/3702.tsv', sep = '\t', header=None, comment="#", skip_blank_lines=True)
    entry_name_df = pd.read_csv('../real_interactions/ara_entry_id_name.csv')
    #print(entry_name_df)
    #print(domain_df)
    data = entry_name_df.drop(['Entry'], axis=1)
    print(data)
    entry_name_dict = data.set_index('From')["Entry Name"].to_dict()
    print(entry_name_dict)
    entry_name_list = []
    for ID in domain_df[0]:
        if ID in entry_name_dict:
            #print(entry_name_dict[ID])
            entry_name_list.append(entry_name_dict[ID])
        else:
            print(ID)
    print(len(entry_name_list))
    print(len(domain_df[0]))
    if len(entry_name_list) == 55713:
        domain_df["Entry Name"] = entry_name_list
    print(domain_df)
    domain_df.to_csv('../real_interactions/Arath_domain_map.tsv', sep='\t', index = False)











if __name__ == '__main__':
    main()



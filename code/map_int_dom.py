
"""
THIS IS a script to map proteins with their protein domains in pfam
and it will sorted by number of proteins in one family
usage: python3 map_int_domain.py
"""
import sys

def readdom(file):
    """function to generate the pfam file into protein: domain_domain

    :param file: str, name of pfam file
    :return: newdict, dictionary {protein: domain_domain}
    """
    f=open(file,'r')
    data=f.readlines()
    f.close()
    dict={}
    for line in data[3:]:
        s = line.split('\t')[:-1]
        prot=s[0]
        dom=s[5]
        # tmp = [domain]
        # dict[prot] = [domain]
        if prot in dict:
            tmp=dict[prot]
        else:
            tmp=[]
        if tmp.count(dom)==0:
            tmp.append(dom)
        tmp.sort()
        dict[prot]=tmp
    newdict={}
    for protein in dict.keys():
        tmp=""
        for domain in dict[protein]:
            tmp=tmp+domain+"_"
        newdict[protein]=tmp[:-1]
    return newdict

def readint(file):
    """function to generate all physical interaction in arabidopsis

    :param file: str, name of file containing all interaction pairs(gene names)
    :return: list of lists, [[a1,b1], [a2,b2],....]
    """
    f=open(file,'r')
    data=f.readlines()
    f.close()
    list=[]
    for line in data:
        s=line.split()
        a=s[0]
        b=s[1]
        if a<b:
            tmp=[a,b]
        else:
            tmp=[b,a]
        if list.count(tmp)==0:
            list.append(tmp)
    return list

prot2dom=readdom(sys.argv[1])
#print(prot2dom)
inter=readint(sys.argv[2])
#print(inter[:10])
countdict={}
# map domains with interaction proteins and count numbers
for i in inter:
    a,b=i
    if a in prot2dom and b in prot2dom:
        dom1=prot2dom[a]
        dom2=prot2dom[b]
        if dom1<dom2:
            tmp=(dom1,dom2)
        else:
            tmp=(dom2,dom1)
        if tmp in countdict:
            c=countdict[tmp]+1
        else:
            c=1
        countdict[tmp]=c
for i in countdict.keys():
    print(i[0],i[1],countdict[i])
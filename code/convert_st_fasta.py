from Bio import SeqIO
from sys import argv
records = SeqIO.parse(argv[1], "stockholm")
count = SeqIO.write(records, argv[2], "fasta")
print("Converted %i records" % count)
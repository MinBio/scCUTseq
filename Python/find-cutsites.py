"""
Get the cutsites in the reference genome (or any other sequence)
"""

from Bio import SeqIO
from Bio.Restriction import *
import re

#specify enzyme
enzyme = NlaIII

#specify reference
ref = "/mnt/AchTeraD/Documents/references/mm10/mm10.fa"

#specify output file
output = "/mnt/AchTeraD/Documents/Projects/scCUTseq/cutsite-distribution/mm10_cutsites.tsv"

#loop through chromosomes
for record in SeqIO.parse(ref, "fasta"):
    sites = enzyme.search(record.seq)
    for x in sites:
        with open(output, "a") as out:
            out.write(re.sub(" .*", "", record.description) +
                "\t" + str(x) + "\t" + str(x + 1) + "\n")
    print(re.sub(" .*", "", record.description) + ": Done")

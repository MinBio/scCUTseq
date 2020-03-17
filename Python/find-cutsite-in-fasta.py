#!/usr/bin/env python3

import sys
from Bio import SeqIO

query = sys.argv[2].lower()

for record in SeqIO.parse(open(sys.argv[1], 'r'), 'fasta'):
    indexes = []
    i = str(record.seq).lower().find(query)
    while i >= 0:
        indexes.append(str(i + 1))
        i = str(record.seq).lower().find(query, i + 1)

    if len(indexes) > 0:
        print(record.id + "\t" + ','.join(indexes))
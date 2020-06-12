#!/usr/bin/env python3

import sys
from Bio import SeqIO

query = sys.argv[2].lower()

for record in SeqIO.parse(open(sys.argv[1], 'r'), 'fasta'):
    indexes = []
    if len(record.seq) > 10000:
	    i = str(record.seq).lower().find(query)
	    while i >= 0:
	        indexes.append(str(i + 1))
	        i = str(record.seq).lower().find(query, i + 1)

	    if len(indexes) > 1:
	    	distances = []
	    	for i in range(len(indexes)-1):
	    		dist = int(indexes[i+1]) - int(indexes[i])
	    		distances.append(dist)
	    	print(record.id + "\t" + ','.join(str(x) for x in distances))
	        #print(record.id + "\t" + ','.join(indexes))

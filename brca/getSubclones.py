import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist


maxdiff = 0.02
method = "single"
minsize = 3
infile = "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/P19254-prostate_integer-cn.tsv"
outfile = "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/P19254-prostate_clones.tsv"

def main(args=None, stdout_file=None):
	# Read in file
	cn = pd.read_csv(infile, delimiter = "\t")

	# Cluster data
	clusters = clustering(cn, maxdiff, method)

	# Select clusters
	clones = selecting(clusters, minsize)

	# Write output
	with open(outfile, "w") as out:
		for i in clones.items():
			out.write("\t".join(str(s) for s in i) + "\n")


def clustering(data, maxdiff, method):
	linkage = hier.linkage(data.T, method=method, metric='hamming', optimal_ordering=True)
	clus = hier.fcluster(linkage, t=maxdiff, criterion='distance')
    return {e : clus[i] for i, e in enumerate(data.columns)}

def selecting(clusters, minsize):
	size = {i : sum(clusters[c] == i for c in clusters) for i in set(clusters.values())}
    return {c : clusters[c] for c in clusters if size[clusters[c]] >= minsize}


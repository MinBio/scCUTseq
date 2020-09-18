import sys
sys.path.append('/home/luukharbers/SCICoNE/pyscicone')
import scicone
import numpy as np

install_path = '/home/luukharbers/SCICoNE/build/'
temporary_outpath = '/mnt/AchTeraD/Documents/Projects/scCUTseq/SCICoNE-testoutput/'
seed = 42 # for reproducibility

sci = scicone.SCICoNE(install_path, temporary_outpath, verbose = True)

# Read in counts
counts = np.loadtxt("/mnt/AchTeraD/Documents/Projects/scCUTseq/hmmcopy/input_matrix/BICRO235_NZ120-500kb-inputmatrix.txt",
    dtype = int)

# Read in chromosome stops as dict
with open("/mnt/AchTeraD/Documents/Projects/scCUTseq/hmmcopy/chromosome-stops-filtered-500kb.tsv") as f:
    chr_stops = dict([line.split() for line in f])

chr_stops = dict((k,int(v)) for k,v in chr_stops.items())


# Normalize counts
normalized_counts = counts / np.sum(np.abs(counts), axis=1).reshape(-1, 1)
normalized_counts *= normalized_counts.shape[1]

# Plot
scicone.plotting.plot_matrix(normalized_counts,
                             cbar_title='Normalized\n  counts', vmax=2,
                             chr_stops_dict=chr_stops,
                             cluster=True)

# Get breakpoints
bps = sci.detect_breakpoints(data=normalized_counts, verbosity=1,
                             window_size=50)

# Plot breakpoints
scicone.plotting.plot_matrix(normalized_counts, bps=bps['segmented_regions'],
                             cbar_title='Normalized\n  counts', vmax=2,
                             chr_stops_dict=chr_stops)

# Infer tree
diploid_tree = sci.learn_tree(normalized_counts, bps['segmented_regions'],
                              cluster_tree_n_iters=10000,
                              cluster=True, n_reps=10, max_tries=4,
                              copy_number_limit=4)


scicone.plotting.plot_matrix(normalized_counts, mode='data',
                             cbar_title='Raw counts', vmax=2,
                             labels=diploid_tree.cell_node_labels)


scicone.plotting.plot_matrix(diploid_tree.outputs['inferred_cnvs'],
                             mode='cnv', cbar_title='CNV', vmax=4,
                             chr_stops_dict=chr_stops,
                             labels=diploid_tree.cell_node_labels,
                             output_path="/home/luukharbers/Desktop/testplot.png")
## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggtree", "tidytree", "treeio")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

tree = read.tree("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/medicc/brca.out/tree_final_A.new")

# Drop diploid and set new root
tree = drop.tip(tree, "diploid")
tree$root.edge = 2

# Plot2
plt = ggtree(tree) + 
  theme_tree() +
  geom_tiplab() +
  geom_nodepoint() +
  geom_rootedge() +
  annotate("text", x = -1.5, y = 3, label = "Diploid")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/subclones-newick_tree",
              height = 6, width = 10)

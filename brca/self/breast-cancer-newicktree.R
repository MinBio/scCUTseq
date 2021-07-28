## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggtree", "tidytree", "treeio")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

tree = read.tree("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/medicc2/brca1/brca1-medicc_dt_final_tree.new")

# Drop diploid and set new root
tree = drop.tip(tree, "N")
tree$root.edge = 2

# Plot2
plt = ggtree(tree) + 
  theme_tree() +
  geom_tiplab() +
  geom_nodepoint() +
  geom_rootedge() +
  annotate("text", x = -2, y = 2.4, label = "Diploid")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/brca1-newicktree",
              height = 6, width = 10)

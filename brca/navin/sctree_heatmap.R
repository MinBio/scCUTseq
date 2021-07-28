## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting ME trees

## Load/install packages
packages = c("data.table", "ape", "ggtree")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")
source("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/plotHeatmap.R")

# Set threads
nthreads = 32

calc_sctree_dists <- function(tree) {
  # getting info for segment annotation
  tbl_tree <- as_tibble(tree)
  
  # maximum length will be the truncal jump
  max_length <-
    max(tbl_tree$branch.length[!is.na(tbl_tree$branch.length)])
  dip_label_height <- which(tbl_tree$branch.length == max_length)
  
  # subsetting a subtree with only nodes after MRCA
  distnodes <- dist.nodes(tree)
  max_mrca_dist <- distnodes %>%
    as_tibble() %>%
    rownames_to_column("node") %>%
    gather(key = "node2", value = "value", -node) %>%
    mutate(node = as.numeric(node),
           node2 = as.numeric(node2)) %>%
    dplyr::filter(node == dip_label_height,
                  node2 >= dip_label_height) %>%
    summarise(max_dist_mrca = max(value)) %>%
    pull(max_dist_mrca)
  
  # return the truncal and branching dist in a df
  distnodes_df <- tibble(truncal_node = dip_label_height,
                         truncal = max_length,
                         branching = max_mrca_dist)
  
  return(distnodes_df)
  
}

plot_sctree <- function(tree,
                        anno_y = 1200,
                        title = NULL) {
  
  dist_nodes <- calc_sctree_dists(tree)
  
  p <-
    ggtree::ggtree(ape::ladderize(tree),
                   ladderize = F,
                   size = .2) +
    geom_tippoint(color = "black", size = 1) +
    theme_tree2(axis.text.x = element_text(size = 15),
                axis.title.x = element_text(size = 15)) +
    # punctuated
    annotate(
      "segment",
      x = 0,
      xend = dist_nodes$truncal - 30,
      y = anno_y,
      yend = anno_y
    ) +
    annotate(
      "text",
      x = dist_nodes$truncal / 2,
      y = anno_y + 30,
      size = 5,
      label = round(dist_nodes$truncal)
    ) +
    annotate(
      "text",
      x = dist_nodes$truncal / 2,
      y = anno_y - 30,
      size = 5,
      label = "truncal"
    ) +
    # MRCA to last node
    annotate(
      "segment",
      x = dist_nodes$truncal + 30,
      xend = dist_nodes$truncal + dist_nodes$branching,
      y = anno_y,
      yend = anno_y
    ) +
    annotate(
      "text",
      x = (dist_nodes$truncal) + (dist_nodes$branching / 2),
      y = anno_y + 30,
      size = 5,
      label = round(dist_nodes$branching)
    ) +
    annotate(
      "text",
      x = (dist_nodes$truncal) + (dist_nodes$branching / 2),
      y = anno_y - 30,
      label = "branching",
      size = 5
    ) +
    annotate(
      "text",
      x = 1400,
      y = dist_nodes$truncal_node,
      label = "diploid",
      size = 5
    ) +
    xlab("manhattan distance")
  
  p
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
  
}
# Load in profiles
profiles = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca1_nondiploid-HQ.rds")

# Make treedata_df and add diploid
tree_data = profiles[, 4:ncol(profiles)]
tree_data[, diploid := 2]
tree_data[, diploid2 := 2]
# run ME tree
tree_me = fastme.bal(dist(t(tree_data), method = "manhattan"))

tree = root.phylo(tree_me,
                  outgroup = "diploid2",
                  resolve.root = T)
tree = drop.tip(tree, tip = c("diploid", "diploid2"))
tree = ladderize(tree)

p = ggtree(tree, ladderize = F, size = .2) +
  geom_tippoint(size = 1, color = "black") +
  theme_tree2(axis.text.x = element_text(size = 15),
              axis.title.x = element_text(size = 15))
# Save sc tree
save_and_plot(p, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/navin/brca2_sctree",
              width = 7, height = 10)

# Save cn heatmap
heatmap = plotHeatmap(profiles[, 4:ncol(profiles)], profiles[, 1:3], dendrogram = T, rasterize = T, linesize = .7)
save_and_plot(heatmap, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/navin/brca1_sc-heatmap",
              width = 24, height = 20)

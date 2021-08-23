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

# Colors
get_colors = function() {
  
  # This function is an unfortunate inheritance from the beginning of this project. Everything uses it
  # so it is was kept.
  
  # subclones
  hues1 = paletteer::paletteer_d("ggsci::default_igv",
                                 n = 35) %>% unclass()
  
  
  colors_vec1 = setNames(hues1, paste0("c", 1:35))
  
  #superclones
  
  colors_g_cl <-
    structure(unclass(paletteer::paletteer_d("yarrr::info", n = 9)),
              names = paste0("s", 1:9))
  colors_g_cl[1] = "#E7A79BFF"
  colors_g_cl[2] = "#90A8C0FF"
  colors_g_cl[6] = "#B4DCF5FF"
  colors_g_cl[7] = "#F2C695FF"
  
  colors_list = list(subclones = colors_vec1,
                     superclones = colors_g_cl)
  
  return(colors_list)
  
}
colors_vector = get_colors()

# Load in profiles
profiles = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca2_nondiploid-HQ.rds")
clones = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/navin/brca2clones.tsv")
cosmic = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/genes/gene-census.tsv")

# Make treedata_df and add diploid
tree_data = melt(profiles, id.vars = c("chr", "start", "end"))
tree_data = merge(tree_data, clones, by.x = "variable", by.y = "id")

# Get consensus profiles
tree_data = tree_data[, .(cn = round(median(value))), by = .(chr, start, end, superclones, subclones)]
tree_data = dcast(tree_data, chr+start+end ~ subclones, value.var = "cn")
tree_data[, diploid := 2]

# run ME tree
tree = fastme.bal(dist(t(tree_data[, 4:ncol(tree_data)]), method = "manhattan"))

tree = root.phylo(tree,
                  outgroup = "diploid",
                  resolve.root = T)
tree = drop.tip(tree, tip = "diploid")
tree = ladderize(tree)

# Order clone info
clones = unique(clones[, 1:2], by = "subclones")
clones[, taxa := subclones]

colors_vector_gg =
  c(colors_vector$subclones, colors_vector$superclones)
p = ggtree(tree, ladderize = F, size = 2) +
  geom_rootedge(rootedge = 300, size = 2)

p = p %<+% clones[, 3:1] +
  geom_tippoint(aes(color = superclones),
                size = 10) +
  geom_tippoint(aes(color = subclones),
                size = 3) +
  scale_color_manual(values = colors_vector_gg) +
  theme(legend.position = "none")

# Save tree
save_and_plot(p, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/navin/brca1_tree", height = 7, width = 9)

# Get tip order
is_tip = tree$edge[, 2] <= length(tree$tip.label)
ordered_tips_index = tree$edge[is_tip, 2]
tree_tips_order = tree$tip.label[ordered_tips_index]

# Plot heatmap
tree_data = tree_data[gtools::mixedorder(tree_data$chr)]
heatmap_dt = tree_data[, 4:(ncol(tree_data)-1), with = F]
combined = plot_grid(p, plotHeatmap(heatmap_dt, tree_data[, 1:3], dendrogram = F, linesize = 11, order = tree_tips_order, rasterize = T), 
                     ncol = 2, rel_widths = c(0.3, 1), axis = "tblr", align = "hv")
save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/navin/brca1_tree_and_heatmap-rast", height = 7, width = 28)

# Get cosmic genes that are amp/del
dt = pivot_longer(tree_data, cols = colnames(tree_data[, 4:ncol(tree_data)]))
setDT(dt)
dt = dt[value != 2, ]
dt[, value := ifelse(value < 2, "Deleted", "Amplified")]

# Setkeys
setkey(dt, chr, start, end)
setkey(cosmic, chr, start, end)

overlap = foverlaps(cosmic, dt)
overlap = overlap[complete.cases(overlap)]
overlap = unique(overlap, by = c("name", "value", "gene"))
overlap = pivot_wider(overlap[, c(1:5, 8)], names_from = "name", values_from = "value")

write.table(overlap, "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/cosmic/brca1_cosmic-overlap.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")

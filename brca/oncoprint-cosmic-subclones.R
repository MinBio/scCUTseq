## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "tidyr", "ComplexHeatmap", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

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

# Get cosmic genes that are amp/del
dt = pivot_longer(tree_data, cols = colnames(tree_data[, 4:ncol(tree_data)]))
setDT(dt)
dt = dt[value != 2, ]
dt[, value := ifelse(value < 2, "Deleted", "Amplified")]

# Setkeys
setkey(dt, chr, start, end)
cosmic[, chr := as.character(chr)]
setkey(cosmic, chr, start, end)

overlap = foverlaps(cosmic, dt)
overlap = overlap[complete.cases(overlap)]
overlap = unique(overlap, by = c("name", "value", "gene"))
overlap = pivot_wider(overlap[, c(1:5, 8)], names_from = "name", values_from = "value")

# Prepare oncoprint matrix
mat = as.matrix(overlap[, 4:ncol(overlap)])
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]

# Remove genes that is the same over all the subclones
select = sapply(1:nrow(mat), function(i) {
  del = sum(mat[i, ] == "Deleted") == ncol(mat)
  amp = sum(mat[i, ] == "Amplified") == ncol(mat)
  return(amp | del)
})
mat = mat[!select, ]

# colors
cols = brewer.pal(3, "Set1")[c(1, 2)]
names(cols) = c("Amplified", "Deleted")

#plot
plt = oncoPrint(mat,
                alter_fun = list(
                 background = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#FFFFFF", col = "black", size = .1),
                 Amplified = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["Amplified"]),
                 Deleted = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["Deleted"])), 
                col = cols, border = "black", show_column_names = T, show_row_names = T, remove_empty_rows = T,
                show_pct = F, row_names_gp = gpar(fontsize = 3))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/oncoprint/brca2-oncoprint",
              height = 12, width = 7)

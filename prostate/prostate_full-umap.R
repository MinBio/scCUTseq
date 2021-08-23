## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "uwot")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load profiles
profiles = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/prostate_p3.rds")
annot = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/P3.tsv", header = F)

samples = data.table(cell = colnames(profiles[, 4:ncol(profiles)]),
                     library = gsub("_.*", "", colnames(profiles[, 4:ncol(profiles)])))
annot = merge(samples, annot, by.x = "library", by.y = "V1", all.x = TRUE)

# Check for equality
stopifnot(all.equal(annot$cell, colnames(profiles[, 4:ncol(profiles)])))


set.seed(45)
total_umap = umap(t(profiles[, 4:ncol(profiles)]), n_neighbors = 24, spread = 1, min_dist = 0, metric = "manhattan")

umap_dt = data.table(x = total_umap[, 1],
                     y = total_umap[, 2],
                     sample = annot$cell,
                     group = annot$V2)

plot = ggplot(umap_dt, aes(x = x, y = y, color = group)) +
  geom_point(size = 1.5) +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "")

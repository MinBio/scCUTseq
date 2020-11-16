## Author: Luuk Harbers
## Date: 2020-11-02
## Script for plotting mapped reads to HS and DM from mixed cell experiments

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/"

rep1_hs = list.files(paste0(base_path, "BICRO220+BICRO217/NZ26+63/hs/"), pattern = "all.tsv$", recursive = T, full.names = T)
rep1_dm = list.files(paste0(base_path, "BICRO220+BICRO217/NZ26+63/dm/"), pattern = "all.tsv$", recursive = T, full.names = T)
rep2_hs = list.files(paste0(base_path, "BICRO217/NZ27/"), pattern = "all_hs.tsv", recursive = T, full.names = T)
rep2_dm = list.files(paste0(base_path, "BICRO217/NZ27/"), pattern = "all_dm.tsv", recursive = T, full.names = T)

sampleinfo = fread(paste0(base_path, "BICRO217/NZ26_NZ27_annot.csv"))

dt = data.table(sample = fread(rep1_hs, select = 1)[[1]],
                rep1_hs = fread(rep1_hs, select = 9)[[1]],
                rep1_dm = fread(rep1_dm, select = 9)[[1]],
                rep2_hs = fread(rep2_hs, select = 9)[[1]],
                rep2_dm = fread(rep2_dm, select = 9)[[1]])

dt = merge(dt, sampleinfo, by.x = "sample", by.y = "BARCODE")
dt[, TYPES := gsub(" POS", "", TYPES)]

# Remove neg ctrls
dt = dt[TYPES != "NEG",]
dt_melt = melt(dt, id.vars = c("sample", "SAMPLE", "TYPES"))
dt_melt[, aligned := ifelse(grepl("hs", variable), "human", "drosophila")]
dt_melt[, replicate := ifelse(grepl("rep1", variable), "Replicate 1", "Replicate 2")]

# Plot values
plot = ggplot(dt_melt, aes(x = aligned, y = value, color = TYPES)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge()) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) + 
  facet_grid(~replicate) +
  scale_color_viridis_d(begin = 0.4) +
  labs(y = "Reads per single cell", x = "Alignment species", color = "Cell type")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/cross-contamination/reads-alignment-replicates",
              width = 7, height = 7)

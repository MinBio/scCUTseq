## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for making normal reference for lrr normalization

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

options(scipen = 999)

binsizes = c(50000, 100000, 175000, 250000, 500000)

lapply(binsizes, function(binsize) {
  rda = paste0("/mnt/AchTeraD/data/BICRO243/MS67/out/dnaobj-psoptim-", binsize, ".Rda")
  rds = paste0("/mnt/AchTeraD/data/BICRO243/MS67/out/dnaobj-", binsize, ".Rds")
  
  # Set QC thresholds
  min_reads = 1e6
  max_spikiness = 0.55
  min_avgreads = 50
  diploid = TRUE
  
  load(rda)
  stats = readRDS(rds)
  bins = stats$binbed[[1]]
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, ..usable_cells]
  
  # Select diploid only
  if(diploid) dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
  
  dnanorm = data.table(dnalrr)
  dnaraw_hq = dnanorm[, colnames(dnanorm) %in% colnames(dt), with = F]
  dnaraw_hq = cbind(bins, dnaraw_hq)
  dnaraw_means = dnaraw_hq[, .(means = rowMeans(.SD)), by = c("chr", "start", "end")]
  
  dnaraw_means[, identifier := factor(paste0(chr, start, end), levels=paste0(chr, start, end))]
  ggplot(dnaraw_means[chr == 8,], aes(x = identifier, y= means)) +
    geom_point() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  # 
  # write.table(dnaraw_means, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/CNA-calling/files/BICRO243_MS67-normal-",
  #                                  binsize, ".tsv"), quote = F, row.names = F, col.names = F, sep = "\t")
})


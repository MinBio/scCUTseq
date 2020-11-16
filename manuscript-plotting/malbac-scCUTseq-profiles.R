## Author: Luuk Harbers
## Date: 2020-10-29
## Script for plotting of selected profiles

## Load/install packages
packages = c("data.table", "tidyverse", "scales", "pbapply", "ggdendro", "GenomicRanges", "ggrastr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")


load("/mnt/AchTeraD/data/BICRO221/out/dnaobj-psoptim-500000.Rda")
data = readRDS("/mnt/AchTeraD/data/BICRO221/out/dnaobj-500000.Rds")

# Select stats
stats = data$stats[[1]]
setDT(stats)

# Select bins
bins = data$binbed[[1]]
setDT(bins)

bins[, bin := seq_along(chr)]
bins[, end_cum := cumsum((end - start) + 1)]
bins[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Make chr_bounds
chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
chr_bounds = chr_bounds %>% 
  mutate(mid = round(min + (max-min) / 2,0),
         end_bp=cumsum(as.numeric(chrlen_bp)), 
         start_bp = end_bp - chrlen_bp, 
         mid_bp = round((chrlen_bp / 2) + start_bp, 0))

#Colors
colors = c("#496bab", "#9fbdd7", "#c1c1c1", "#e9c47e", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# Plot profile plots
invisible(lapply(samples, function(sample) {
  dt = cbind(bins, psoptim[, sample], 2^dnalrr[,sample]*psoptim.par[, sample])
  setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "cn", "raw"))
  
  dt[, col := ifelse(cn < 11, as.character(cn), "10+")]
  dt[, col := factor(col, levels = c(as.character(0:10), "10+"))]
  
  # Sample stats
  stat_string = paste0(sample, 
                       " | reads: ", stats[cell == sample,]$reads,
                       " | avg reads/bin: ", as.integer(stats[cell == sample]$mean),
                       " | spikiness: ", round(stats[cell == sample,]$spikiness, 2))
  
  # save plot
  plt = ggplot(dt, aes(x = bin)) +
    geom_point(aes(y = raw, color = col), size = 0.7) +
    geom_point(aes(y = cn), size = 1) +
    scale_color_manual(values = colors, drop = F) + 
    scale_y_continuous(labels=comma_format(accuracy = 1), breaks = pretty_breaks(6)) + 
    scale_x_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid) + 
    labs(y = "Copy Number", x = "", subtitle = stat_string, color = "") + 
    geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) + 
    theme(legend.position = "right",
          axis.ticks.x = element_blank())
  # Save plot
  save_and_plot(plt, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac-scCUTseq/BICRO221/", sample), 
         width = 14, height = 5)
}))

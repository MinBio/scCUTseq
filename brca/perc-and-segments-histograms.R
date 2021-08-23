## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for saving rds of HQ profiles

## Load/install packages
packages = c("data.table", "pbapply", "tidyr", "GenomicRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
nthreads = 32

cn = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca2.rds")

total_size = sum(cn$end - cn$start)
cn_m = melt(cn, id.vars = c("chr", "start", "end"))

cn_m[, length := end-start]
amp_perc = cn_m[value > 2, .(AMP = sum(length) / total_size), by = variable]
del_perc = cn_m[value < 2, .(DEL = sum(length) / total_size), by = variable]

perc = merge(amp_perc, del_perc)
# Add missing cells
perc = rbind(perc, 
               data.table(variable = unique(cn_m$variable)[!unique(cn_m$variable) %in% perc$variable],
                          AMP = 0,
                          DEL = 0))

perc = melt(perc, id.vars = "variable", variable.name = "alteration")

# Plot
plt1 = ggplot(perc, aes(x = value, fill = alteration)) +
  geom_histogram() +
  facet_wrap(~alteration) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Number of Cells", x = "Percentage of genome altered", fill = "")

plt2 = ggplot(perc, aes(x = value, fill = alteration)) +
  geom_histogram(alpha = .55, position = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Number of Cells", x = "Percentage of genome altered", fill = "")

# Get number of segments
cn_m[, start_alt := 1:.N, by = .(chr, variable)]
cn_m[, end_alt := 1:.N+1, by = .(chr, variable)]

segments = cn_m[, as.data.table(reduce(IRanges(start_alt, end_alt))), by = .(chr, variable, value)]
amp_seg = segments[value > 2, .(AMP = .N), by = variable]
del_seg = segments[value < 2, .(DEL = .N), by = variable]
counts = merge(amp_seg, del_seg)

# Add missing cells
counts = rbind(counts, 
               data.table(variable = unique(cn_m$variable)[!unique(cn_m$variable) %in% counts$variable],
               AMP = 0,
               DEL = 0))

counts_m = melt(counts)
setnames(counts_m, c("variable", "alteration", "value"))

plt3 = ggplot(counts_m, aes(x = value, fill = alteration)) +
  geom_histogram() +
  facet_wrap(~alteration) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Number of Cells", x = "Number of segments altered", fill = "")

plt4 = ggplot(counts_m, aes(x = value, fill = alteration)) +
  geom_histogram(alpha = .55, position = "identity") +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Number of Cells", x = "Number of segments altered", fill = "")

# Save plots
save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/histograms/brca2_percentage-altered-split",
              width = 12, height = 6)
save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/histograms/brca2_percentage-altered",
              width = 6, height = 6)
save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/histograms/brca2_num-segments-split",
              width = 12, height = 6)
save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/histograms/brca2_num-segments",
              width = 6, height = 6)

## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 


## Load/install packages
packages = c("data.table", "ggplot2", "GenomicRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/BICRO253/"
nthreads = 40

rda = list.files(base_path, pattern = "-500000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-500000.Rds", recursive = T, full.names = T)

# Set thresholds for HQ profiles
min_reads = 3e5
max_spikiness = 0.55
min_avgreads = 50

rds = readRDS(rds)
load(rda)
stats = rds$stats[[1]]
setDT(stats)
dt = data.table(psoptim)

usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
dt = dt[, ..usable_cells]

# Select diploid only
dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
bins = rds$binbed[[1]]
total = cbind(bins, dt)
total = total[chr != "chrY",]

# Get segments
segments = lapply(colnames(total)[4:ncol(total)], function(x) {
  res = total[, c("chr", "start", "end", x), with = F]
  res = res[res[[x]] < 2 | res[[x]] > 2,]
  if(nrow(res > 0)) {
    gr = makeGRangesFromDataFrame(res, keep.extra.columns = T)
    red = reduce(gr)
    dt = as.data.table(red)
    tot = data.table(sample = x, start = dt$start, end = dt$end, cn = res[start %in% dt$start,][[4]])
    return(tot)
  }
})
total = rbindlist(segments)
total[, alteration := ifelse(cn < 2, "Deletion", "Amplification")]
add = colnames(dt)
add = data.table(sample = rep(add, 2), alteration = c(rep("Amplification", length(add)), rep("Deletion", length(add))), N = 0)

counts = total[, .N, by = .(sample, alteration)]
counts = merge(add, counts, by = c("sample", "alteration"), all.x = T)
counts[is.na(`N.y`), `N.y` := 0]

number = counts[, .N, by = .(alteration)]
number_1 = counts[`N.y` > 5, .N, by = .(alteration)]

plt = ggplot(counts, aes(x = alteration, y = `N.y`, color = alteration)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  geom_text(data=number, aes(y = -2, label = paste0("n = ", N)), hjust = 0.5, color = "black") +
  scale_color_brewer(palette = "Set1") +
  labs(y = "Number of alterations", x = "") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/mouse/amp-del-segments-boxplot",
              height = 7, width = 7)

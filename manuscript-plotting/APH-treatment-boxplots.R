## Author: Luuk Harbers
## Date: 2020-10-29
## Script for plotting boxplots of CNA burden after APH treatment

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/BICRO230+232+233/"

rda = list.files(base_path, pattern = ".Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = ".Rds", recursive = T, full.names = T)
libraries = basename(list.dirs(base_path, recursive = F))


data = lapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  return(list(data.table(psoptim), stats))
})

# Set thresholds
min_reads = 3e5
max_spikiness = 0.7
min_avgreads = 50

res = lapply(1:length(data), function(lib) {
  # Get usable cells and filter dt
  usable_cells = data[[lib]][[2]][reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = data[[lib]][[1]][, usable_cells, with = F]
  return(dt)
})

cna = lapply(res, function(lib) {
  amp = apply(lib, 2, function(x) sum(x > 2)) / nrow(lib) * 100
  del = apply(lib, 2, function(x) sum(x < 2)) / nrow(lib) * 100
  
  return(data.table(sample = names(amp), amp = amp, del = del))
})

# Plot
cna[[1]][, treatment := "NT"]
cna[[2]][, treatment := "72hr APH"]

cna_total = rbindlist(cna)
cna_total[, CNA := amp + del]
cna_total = melt(cna_total[, .(sample, treatment, CNA)], id.vars = c("sample", "treatment"))
cna_total[, treatment := factor(treatment, levels = c("NT", "72hr APH"))]

obs = cna_total[, .N, by = treatment]


plt = ggplot(cna_total[variable == "CNA"], aes(x = treatment, y = value, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) + 
  stat_compare_means(method = "wilcox.test", paired = F, comparisons = list(c("72hr APH", "NT"))) +
  geom_text(data = obs, aes(y = 0, label = paste("n =", N)), color = "black") +
  labs(y = "Genome altered (%)", x = "Treatment") + 
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme(legend.position = "none")
  
save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/APH-treatment/NZ75-NZ76_72hrAPH-treatment",
              width = 7, height = 7)


## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting pearson's correlation for different malbac scalings  

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/BICRO226+227/"

ref_rda = "/mnt/AchTeraD/data/BICRO188/XZ244/out/dnaobj-psoptim-500000.Rda"
rda = list.files(base_path, pattern = ".Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = ".Rds", recursive = T, full.names = T)
libraries = basename(list.dirs(base_path, recursive = F))
libraries_match = c("1:50", "1:100", "1:200", "1:500")

getref = function(ref_rda) {
  load(ref_rda)
  return(data.table(psoptim[, 2]))
}
ref_profile = getref(ref_rda)

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
min_avgreads = 30

res = lapply(1:length(data), function(lib) {
  # Get usable cells and filter dt
  usable_cells = data[[lib]][[2]][reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = data[[lib]][[1]][, usable_cells, with = F]
  
  # Get pairwise pearson's correlation
  cors = apply(dt, 2, function(x) {
    cor(x, ref_profile)  
  })
  # pairwise = cor(dt)
  # pairwise[upper.tri(pairwise, diag = T)] = NA
  # pairwise = reshape2::melt(pairwise, na.rm = T)
  # pairwise$sample = libraries_match[lib]
  return(data.table(cor = cors, sample = libraries_match[lib]))
})

# bind results
total = rbindlist(res)

# Set factor levels
total[, sample := factor(sample, levels = libraries_match)]
comparisons = list(c("1:50", "1:100"),
                   c("1:100", "1:200"),
                   c("1:200", "1:500"))

obs = total[, .N, by = sample]
# Plot
plot = ggplot(total, aes(x = sample, y = cor)) +
  geom_violin(aes(fill = sample, color = sample)) +
  geom_boxplot(width = 0.05, outlier.size = 1) +
  stat_compare_means(comparisons = comparisons, paired = F,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  geom_text(data = obs, aes(y = 0, label = paste("n =", N))) +
  scale_fill_viridis_d(begin = 0.4) +
  scale_color_viridis_d(begin = 0.4) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
  labs(y = "Pearson's Correlation to bulk SKBR3 profile", x = "MALBAC scaling") + 
  theme(legend.position = "none") 

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac_scaling/pearson-to-bulk-BICRO226_227-malbac_scaling",
              height = 7, width = 7)

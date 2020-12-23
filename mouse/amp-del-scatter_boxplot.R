## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 


## Load/install packages
packages = c("data.table", "ggplot2", "pbapply", "ggdendro")
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


# Loop through libraries
data = pblapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, ..usable_cells]
  
  # Select diploid only
  dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
  # Return dt
  return(dt)
}, cl = nthreads)

# Prepare for plotting heatmap
stats = readRDS(rds[1])
bins = stats$binbed[[1]]
setDT(bins)

ind = which(bins$chr != "chrY")
bins = bins[ind,]
total = do.call(cbind, data)
dt = data.table(cbind(bins, total[ind,]))
dt = dt[chr != "chrY",]

means = dt[, .(mean = rowMeans(.SD)), by = .(chr, start, end)]
blacklist = c(3790:3839, 1720:1779, 3016:3026)

dt[blacklist, ] = NA
dt = dt[complete.cases(dt)]

cnas = lapply(dt[, 4:ncol(dt)], function(x) {
  dels = sum(x < 2)
  amps = sum(x > 2)
  return(data.table(amp = amps/nrow(dt), del = dels/nrow(dt)))
})

cnas = rbindlist(cnas)
cnas[, sample := colnames(dt)[4:ncol(dt)]]
setnames(cnas, c("Amplification", "Deletion", "sample"))
cnas_m = melt(cnas, by = "sample")

plt1 = ggplot(cnas, aes(x = Amplification, y = Deletion)) +
  geom_point() +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(y = "Deletion (%)", x = "Amplification (%)")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/mouse/amp-del-scatterplot",
              height = 7, width = 7)

number = cnas_m[, .N, by = .(variable)]
number_1 = cnas_m[value > 0.01, .N, by = .(variable)]
plt2 = ggplot(cnas_m, aes(x = variable, y = value, color = variable)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2) +
  geom_text(data=number, aes(y = -0.009, label = paste0("n = ", N)), hjust = 0.5, color = "black") +
  geom_text(data=number_1, aes(y = 0.1, label = paste0("n (>1%) = ", N)), hjust = 0.5, color = "black") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_brewer(palette = "Set1") +
  labs(y = "Percentage altered", x = "") +
  theme(legend.position = "none")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/mouse/amp-del-boxplot",
              height = 7, width = 7)

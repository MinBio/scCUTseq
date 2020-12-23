## Author: Luuk Harbers
## Date: 2020-10-30
## Script for plotting correlation heatmap

## Load/install packages
packages = c("data.table", "ggplot2", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
rda = list.files("/mnt/AchTeraD/data/", pattern = "-500000.Rda", recursive = T, full.names = T)
rda = rda[grepl("BICRO221/NEB|NZ39|NZ40", rda)]
rds = list.files("/mnt/AchTeraD/data/", pattern = "-500000.Rds", recursive = T, full.names = T)
rds = rds[grepl("BICRO221/NEB|NZ39|NZ40", rds)]

libs = c("scCUTseq-live", "scCUTseq-fixed", "sMALBAC")

# Set QC thresholds
min_reads = 1e6
max_spikiness = 0.55
min_avgreads = 50

total = pblapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, ..usable_cells]

  # Set colnames to include library
  setnames(dt, paste0(libs[lib], "_", colnames(dt)))
  if(lib > 2) {
    setnames(dt, paste0(libs[lib], c("-fixed (cell 3)", "-live (cell 1)", "-fixed (cell 1)", "-live (cell 2)", 
                                     "-live (cell 3)", "-fixed (cell 2)")))
  }
  # Return dt
  return(dt)
})

dt = do.call(cbind, total)

#dt = dt[, grepl("MS35|MS36|MS37|MS38|NZ47|NZ49|NZ50", colnames(dt)), with = F]

res = cor(dt)

res_m = melt(res, na.rm = T)
setDT(res_m)

res_m = res_m[grepl("scCUTseq", Var1) & grepl("sMALBAC", Var2), ]
live = res_m[grepl("live", Var1) & grepl("live", Var2), ]
fixed = res_m[grepl("fixed", Var1) & grepl("fixed", Var2), ]

plt1 = ggplot(live, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_y_discrete("") + 
  scale_x_discrete("", position = "top") +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_viridis("Pearson's\ncorrelation", option="B", begin = 0.75, direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
        axis.line = element_blank())


plt2 = ggplot(fixed, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_y_discrete("") + 
  scale_x_discrete("", position = "top") +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_viridis("Pearson's\ncorrelation", option="B", begin = 0.75, direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
        axis.line = element_blank())

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac-scCUTseq/sMALBAC-vs-scCUTseq-correlation_heatmap-live",
              height = 7, width = 7.5)
save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac-scCUTseq/sMALBAC-vs-scCUTseq-correlation_heatmap-fixed",
              height = 5, width = 10)

## Author: Luuk Harbers
## Date: 2020-10-30
## Script for plotting Lorenz curves

## Load/install packages
packages = c("data.table", "ineq")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
files = list.files("/mnt/AchTeraD/data/BICRO229/NEB/bincounts", full.names = T)
files = c(files, list.files("/mnt/AchTeraD/data/BICRO221/bincounts/", full.names = T))
files = files[grepl("MS35|MS36|MS37|MS38|NZ47|NZ49|NZ50", files)]

data = lapply(files, function(file) {
  dt = fread(file)
  # Get lorenz curve points
  lc = Lc(dt[[1]])
  
  return(data.table(l = lc$L, p = lc$p, sample = basename(file)))
})

plots = lapply(data, function(x) {
  x[, sample := gsub("_S.*|_500000.*", "", sample)]
  x[, type := ifelse(grepl("NZ", sample), "NEBNext 1:200", "NEBNext 1:1")]
  
  # Plot lorenz curve
  plot = ggplot(x, aes(x = p, y = l)) +
    geom_abline(slope = 1, linetype = 2) +
    geom_path(aes(group = sample), size = 1.5, color = "red") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = paste0(x$sample[1], " - ", x$type[1]), 
         y = "Cumulative fraction of total reads", 
         x = "Cumulative fraction of genome", color = "") +
    coord_equal()
  
  save_and_plot(plot, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/lorenz-curve/NEBNext-malbac_scaling-",
                             x$sample[1], "_", x$type[1]),
                height = 7, width = 7)
})


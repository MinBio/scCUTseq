## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ineq", "colorblindr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

raw = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/lorenz-counts-500kb.rds")

# Give names
setnames(raw, c("NEBNext - live 1:1", "NEBNext - fixed 1:1", "NEBNext - live 1:200", 
                "NEBNext - fixed 1:200", "scCUTseq - Cell 1", "scCUTseq - Cell 2", "Bulk CUTseq"))

lorenz = lapply(colnames(raw), function(sample) {
  # Get lorenz curve points
  lc = Lc(raw[[sample]])
  
  return(data.table(l = lc$L, p = lc$p, sample = sample))
})

lorenz = rbindlist(lorenz)

plot = ggplot(lorenz, aes(x = p, y = l, color = sample)) +
  geom_abline(slope = 1, size = 1.25) +
  geom_path(aes(group = sample), size = 1.25) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_npg() +
  labs(y = "Cumulative fraction of total reads", 
       x = "Cumulative fraction of genome", color = "") +
  coord_equal()

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/technical-validation/lorenz/lorenz-comparison",
              width = 8, height = 6)

# Get Gini coefficient
gini = apply(raw, 2, function(x) ineq(x, type = "Gini"))


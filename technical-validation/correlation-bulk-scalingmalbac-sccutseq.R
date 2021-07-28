## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

raw = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/lorenz-counts-500kb.rds")

# Give names
setnames(raw, c("NEBNext - live 1:1", "NEBNext - fixed 1:1", "NEBNext - live 1:200", 
                "NEBNext - fixed 1:200", "scCUTseq - Cell 1", "scCUTseq - Cell 2", "Bulk CUTseq"))

res = cor(raw)
dt = data.table(V1 = rownames(res)[row(res)], 
                V2 = colnames(res)[col(res)], 
                Correlation = c(res))

dt[, V1 := factor(V1, levels = rev(c("Bulk CUTseq", "NEBNext - live 1:1", "NEBNext - fixed 1:1", "NEBNext - live 1:200", 
                                 "NEBNext - fixed 1:200", "scCUTseq - Cell 2", "scCUTseq - Cell 1")))]
dt[, V2 := factor(V2, levels = c("Bulk CUTseq", "NEBNext - live 1:1", "NEBNext - fixed 1:1", "NEBNext - live 1:200", 
                                 "NEBNext - fixed 1:200", "scCUTseq - Cell 2", "scCUTseq - Cell 1"))]
plt = ggplot(dt, aes(x = V1, y = V2, fill = Correlation)) +
  geom_tile() +
  scale_y_discrete("") + 
  scale_x_discrete("", position = "top") +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_viridis("Pearson's\ncorrelation", option="B", begin = 0.75, direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
        axis.line = element_blank())

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/technical-validation/correlation/bulk-malbacscaling-scCUTseq-cor",
              width = 7, height = 6)

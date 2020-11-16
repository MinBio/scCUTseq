## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting duplication rate for different malbac scalings  

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/BICRO226+227/"

mapped = list.files(base_path, pattern = "_all.tsv", recursive = T, full.names = T)
dedup = list.files(base_path, pattern = "_dedup.tsv", recursive = T, full.names = T)
mapped = c("/mnt/AchTeraD/data/BICRO217/NZ28/all_hs.tsv", mapped)
dedup = c("/mnt/AchTeraD/data/BICRO217/NZ28/dedup_hs.tsv", dedup)

libraries = data.table(library = basename(c("NZ28", list.dirs(base_path, recursive = F))), 
                       scaling = c("CUTseq", "1:50", "1:100", "1:200", "1:500"))

dup_rate = lapply(1:length(mapped), function(lib) {
  map = fread(mapped[lib])
  dup = fread(dedup[lib])
  
  dt = data.table(map = map$`#Mapped`, dup = dup$`#Mapped`)
  dt[, rate := (1-(dup/map)) * 100]
  dt[, scaling := libraries$scaling[lib]]
  return(dt[dup > 50e3,])
})

dup_rate = rbindlist(dup_rate)
dup_rate[, scaling := factor(scaling, levels = libraries$scaling)]

comparisons = list(c("1:50", "1:100"),
                   c("1:100", "1:200"),
                   c("1:200", "1:500"))

obs = dup_rate[, .N, by = scaling]

plot1 = ggplot(dup_rate, aes(x = scaling, y = rate)) +
  geom_violin(aes(fill = scaling, color = scaling)) +
  geom_boxplot(width = 0.05, outlier.size = 1) +
  geom_text(data = obs, aes(y = 10, label = paste("n =", N))) +
  scale_fill_viridis_d(begin = 0.4) +
  scale_color_viridis_d(begin = 0.4) +
  labs(y = "Duplication rate (%)", x = "MALBAC scaling") +
  theme(legend.position = "none")
  
save_and_plot(plot1, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac_scaling/duplication_rate-CUTseq-malbac_scaling",
              height = 7, width = 7)

plot2 = ggplot(dup_rate, aes(x = scaling, y = rate)) +
  geom_jitter(width = 0.15, color = "#8b0000") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, size = .75) + 
  geom_text(data = obs, aes(y = 10, label = paste("n =", N)), color = "black") +
  scale_fill_viridis_d(begin = 0.4) +
  scale_color_viridis_d(begin = 0.4) +
  labs(y = "Duplication rate (%)", x = "MALBAC scaling") +
  theme(legend.position = "none")

save_and_plot(plot2, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac_scaling/duplication_rate-CUTseq-malbac_scaling_dotplot",
              height = 7, width = 7)

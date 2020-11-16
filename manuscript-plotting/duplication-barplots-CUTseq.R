## Author: Luuk Harbers
## Date: 2020-10-29
## Script for plotting duplication rates of one library (old preprocessing)

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/BICRO217/NZ28/"

all = fread(paste0(base_path, "all_hs.tsv"))
dedup = fread(paste0(base_path, "dedup_hs.tsv"))

# Combine
dt = data.table(sample = all$Sample, mapped = all$`#Mapped`, deduplicated = dedup$`#Mapped`)
dt[, duplication_rate := (1 - deduplicated / mapped) * 100]
setorder(dt, -mapped)

# Melt and set factors
dt_melt = melt(dt[mapped > 50e3, ], id.vars = "sample")
dt_melt[, sample := factor(sample, levels = dt$sample)]

# Plot
plt = ggplot(dt_melt[variable != "duplication_rate"], aes(x = sample, y = value, fill = variable)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Reads", x = "Cells (>50K reads)", fill = "") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/duplication-plots/CUTseq-on-singlecells-NZ28-duplication-barplot",
              width = 8, height = 5)

obs = data.table(N = dt[mapped > 50e3, .N])

# violin plot
plot = ggplot(dt[mapped > 30e3,], aes(x = "CUTseq on single cells (>50K reads)", y = duplication_rate)) +
  geom_violin() +
  geom_boxplot(width = 0.05, outlier.size = 1.5, size = .7) +
  geom_text(data = obs, aes(y = 70, label = paste("n =", N))) +
  scale_y_continuous(limits = c(70, 82)) +
  labs(y = "Duplication rate (%)", x = "") +
  theme(legend.position = "none")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/duplication-plots/CUTseq-on-singlecells-NZ28-duplication-violinplot",
              width = 5, height = 7)

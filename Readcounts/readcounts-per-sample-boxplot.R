#generate read counts plots scCUTseq

require(data.table)
require(RColorBrewer)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

all <- fread("/mnt/AchTeraD/data/BICRO213/all.tsv", select = c(1, 9))
dedup <- fread("/mnt/AchTeraD/data/BICRO213/dedup.tsv", select = c(1, 9))
setnames(all, "#Mapped", "mapped")
setnames(dedup, "#Mapped", "deduplicated")

data <- merge(all, dedup)

data[, mapped := mapped / 1e6]
data[, deduplicated := deduplicated / 1e6]
data <- melt(data)

plt <- ggplot(data, aes(x = variable, y = value, color = variable)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = c(0.3, 0.05), linetype = 2) +
  scale_color_brewer(palette="Set1") +
  labs(y = "Reads (millions)",
       x = "",
       color = "")


save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/read-distribution/BICRO213_read-distribution", height = 7, width = 7)

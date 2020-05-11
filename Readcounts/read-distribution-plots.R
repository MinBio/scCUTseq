require(data.table)
require(RColorBrewer)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#load in all files
annot <- fread("/mnt/AchTeraD/data/BICRO217/NZ26_NZ27_annot.csv", select = 2:3)
setnames(annot, c("Sample", "cell"))

# MS20_all <- fread("/mnt/AchTeraD/data/BICRO210/MS20/MS20.all.tsv", select = c(1, 9))
# MS20_dedup <- fread("/mnt/AchTeraD/data/BICRO210/MS20/MS20.dedup.tsv", select = c(1, 9))

all_hs <- fread("/mnt/AchTeraD/data/BICRO217/NZ27/all_hs.tsv", select = c(1, 9))
dedup_hs <- fread("/mnt/AchTeraD/data/BICRO217/NZ27/dedup_hs.tsv", select = c(1, 9))
all_s2 <- fread("/mnt/AchTeraD/data/BICRO217/NZ27/all_dm.tsv", select = c(1, 9))
dedup_s2 <- fread("/mnt/AchTeraD/data/BICRO217/NZ27/dedup_dm.tsv", select = c(1, 9))

setnames(all_hs, c("Sample", "all_hs"))
setnames(dedup_hs, c("Sample", "dedup_hs"))
setnames(all_s2, c("Sample", "all_s2"))
setnames(dedup_s2, c("Sample", "dedup_s2"))

#merge

all <- lapply(ls(pattern = "all|dedup"), get)
all <- lapply(all, function(x) merge(x, annot, all=TRUE))

merged <- Reduce(function(...) merge(..., all=TRUE), all)
merged <- merged[, c(1, 2, 4, 6, 8, 9)]
setnames(merged, "cell.y", "cell")
merged <- melt(merged, id.vars = c("Sample", "cell"))

merged$value <- merged$value / 1e6
merged[, cell := gsub(" POS", "", cell)]
merged[, variable := factor(variable, levels = c("all_hs", "dedup_hs", "all_s2", "dedup_s2"))]

plt <- ggplot(merged, aes(x=variable, y=value, color = cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .2)) +
  geom_hline(yintercept = 0.3, linetype = 2) +
  labs(y = "reads (millions)",
       x = "",
       color = "") +
  scale_color_manual(values = brewer.pal(4, "Set1")) +
  theme(legend.position = "top")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/read-distribution/BICRO217/NZ27_all-cells", height = 8, width = 8)
  
require(data.table)
require(RColorBrewer)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#load in all files
annot <- fread("/mnt/AchTeraD/data/BICRO210/annotation-MS22-MS24.csv")
setnames(annot, c("cell", "Sample"))
annot[Sample == "CGGATCGA"]$cell <- "SKBR3"

# MS20_all <- fread("/mnt/AchTeraD/data/BICRO210/MS20/MS20.all.tsv", select = c(1, 9))
# MS20_dedup <- fread("/mnt/AchTeraD/data/BICRO210/MS20/MS20.dedup.tsv", select = c(1, 9))

MS22_all <- fread("/mnt/AchTeraD/data/BICRO210/MS22/all.tsv", select = c(1, 9))
MS22_dedup <- fread("/mnt/AchTeraD/data/BICRO210/MS22/dedup.tsv", select = c(1, 9))
# MS22_all_s2 <- fread("/mnt/AchTeraD/data/BICRO210/MS22/s2/MS22_s2_all.tsv", select = c(1, 9))
# MS22_dedup_s2 <- fread("/mnt/AchTeraD/data/BICRO210/MS22/s2/MS22_s2_dedup.tsv", select = c(1, 9))

MS24_all <- fread("/mnt/AchTeraD/data/BICRO210/MS24/MS24_all.tsv", select = c(1, 9))
MS24_dedup <- fread("/mnt/AchTeraD/data/BICRO210/MS24/MS24_dedup.tsv", select = c(1, 9))
# MS24_all_s2 <- fread("/mnt/AchTeraD/data/BICRO210/MS24/s2/MS24_s2_all.tsv", select = c(1, 9))
# MS24_dedup_s2 <- fread("/mnt/AchTeraD/data/BICRO210/MS24/s2/MS24_s2_dedup.tsv", select = c(1, 9))


# setnames(MS20_all, c("Sample", "MS20_all"))
# setnames(MS20_dedup, c("Sample", "MS20_dedup"))

setnames(MS22_all, c("Sample", "MS22_all"))
setnames(MS22_dedup, c("Sample", "MS22_dedup"))
# setnames(MS22_all_s2, c("Sample", "MS22_s2_all"))
# setnames(MS22_dedup_s2, c("Sample", "MS22_s2_dedup"))
setnames(MS24_all, c("Sample", "MS24_all"))
setnames(MS24_dedup, c("Sample", "MS24_dedup"))
# setnames(MS24_all_s2, c("Sample", "MS24_s2_all"))
# setnames(MS24_dedup_s2, c("Sample", "MS24_s2_dedup"))


all <- lapply(ls(pattern = "MS"), get)
all <- lapply(all, function(x) merge(x, annot, all=TRUE))

merged <- Reduce(function(...) merge(..., all=TRUE), all)
merged <- merged[, c(1, 2, 4, 6, 8, 9)]
setnames(merged, "cell.y", "cell")
merged <- melt(merged, id.vars = c("Sample", "cell"))

merged$value <- merged$value / 1e6

plt <- ggplot(merged, aes(x=variable, y=value, color = cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .2)) +
  geom_hline(yintercept = 0.3, linetype = 2) +
  labs(y = "reads (millions)",
       x = "",
       color = "") +
  scale_color_manual(values = brewer.pal(4, "Set1")) +
  theme(legend.position = "top")

#save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/read-distribution/BICRO210-MS22+24_human-read-distribution", height = 8, width = 8)

merged$value <- merged$value * 1e6
duplication <- dcast(merged, Sample + cell ~ variable)
duplication$MS22_percent <- (duplication$MS22_all - duplication$MS22_dedup) / duplication$MS22_all * 100
duplication$MS24_percent <- (duplication$MS24_all - duplication$MS24_dedup) / duplication$MS24_all * 100
dup <- melt(duplication[, c(1:2, 7:8)], id.vars = c("Sample", "cell"))

plt <- ggplot(dup[cell != "S2" & cell != "neg"], aes(x=variable, y=value, color = cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .2)) +
  labs(y = "Duplicated (%)",
       x = "",
       color = "") +
  scale_color_manual(values = brewer.pal(4, "Set1")) +
  theme(legend.position = "top")

#save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/read-distribution/BICRO210-MS22+24_human-duplication", height = 8, width = 8)

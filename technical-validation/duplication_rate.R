## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for calculating duplication rate

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

duplicate_files = data.table(sample = c("CUTseq", "1:50", "1:100", "1:200", "1:500"),
                             total = c("/mnt/AchTeraD/data/BICRO217/NZ28/NZ28_S3_all.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ72+76/all.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ73+77/all.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ74+78/all.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ75+79/all.tsv"),
                             dedup = c("/mnt/AchTeraD/data/BICRO217/NZ28/NZ28_S3_dedup.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ72+76/dedup.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ73+77/dedup.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ74+78/dedup.tsv",
                                       "/mnt/AchTeraD/data/BICRO226+227/NZ75+79/dedup.tsv"))

res = lapply(1:nrow(duplicate_files), function(i) {
  total = fread(duplicate_files[i, total])
  dedup = fread(duplicate_files[i, dedup])
  
  duplication_percentage = 1 - (dedup[total$`#Mapped` > 50000, `#Mapped`] / total[`#Mapped` > 50000, `#Mapped`])

  return(data.table(library = duplicate_files[i, sample],
                    duplication = duplication_percentage))
})
res = rbindlist(res)
res[, library := factor(library, levels = duplicate_files$sample)]
counts = res[, .N, by = library]

plt = ggplot(res, aes(x = library, y = duplication)) +
  geom_violin(aes(fill = library, color = library)) + 
  geom_boxplot(width = .075) +
  scale_fill_viridis_d(begin = .3) +
  scale_color_viridis_d(begin = .3) +
  scale_y_continuous(labels = scales::label_percent()) +
  geom_text(data = counts, aes(y = 0.1, x = library, label = paste0("n = ", N))) +
  labs(y = "Duplication rate", x = "MALBAC scaling") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/technical-validation/duplication/CUTseq-vs-malbac-scaling",
              width = 7, height = 7)

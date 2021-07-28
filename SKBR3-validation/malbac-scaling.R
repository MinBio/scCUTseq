## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for analysis of scale down malbac scCUTse

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

ref = readRDS("/mnt/AchTeraD/data/BICRO188/XZ244/cnv/500000/cnv.rds")
mb50 = readRDS("/mnt/AchTeraD/data/BICRO226+227/NZ72+76/cnv/500000/cnv.rds")
mb100 = readRDS("/mnt/AchTeraD/data/BICRO226+227/NZ73+77/cnv/500000/cnv.rds")
mb200 = readRDS("/mnt/AchTeraD/data/BICRO226+227/NZ74+78/cnv/500000/cnv.rds")
mb500 = readRDS("/mnt/AchTeraD/data/BICRO226+227/NZ75+79/cnv/500000/cnv.rds")

# Select HQ cells
ref_cn = ref$copynumber$TGATGCGC
mb50_cn = mb50$copynumber[, mb50$stats[classifier_prediction == "good", sample], with = F]
mb100_cn = mb100$copynumber[, mb100$stats[classifier_prediction == "good", sample], with = F]
mb200_cn = mb200$copynumber[, mb200$stats[classifier_prediction == "good", sample], with = F]
mb500_cn = mb500$copynumber[, mb500$stats[classifier_prediction == "good", sample], with = F]
mb50_cn = mb50$copynumber
mb100_cn = mb100$copynumber
mb200_cn = mb200$copynumber
mb500_cn = mb500$copynumber

# Comparison against ref
mb50_cor = data.table(sample = "1:50", pearson = as.vector(cor(mb50_cn, ref_cn)))
mb100_cor = data.table(sample = "1:100", pearson = as.vector(cor(mb100_cn, ref_cn)))
mb200_cor = data.table(sample = "1:200", pearson = as.vector(cor(mb200_cn, ref_cn)))
mb500_cor = data.table(sample = "1:500", pearson = as.vector(cor(mb500_cn, ref_cn)))

res = rbindlist(list(mb50_cor, mb100_cor, mb200_cor, mb500_cor))
res[, sample := factor(sample, levels = c("1:50", "1:100", "1:200", "1:500"))]
obs = res[, .N, by = sample]

# Plot
plt1 = ggplot(res, aes(x = sample, y = pearson)) +
  geom_violin(aes(fill = sample, color = sample)) +
  geom_boxplot(width = .1) +
  geom_text(data = obs, aes(y = 0, label = paste("n =", N))) +
  scale_fill_viridis_d(begin = .4) +
  scale_color_viridis_d(begin = .4) +
  scale_y_continuous(limit = c(0, 1.3), breaks = seq(0, 1, .2), labels = seq(0, 1, .2)) +
  labs(y = "Pearson's Correlation to bulk CUTseq", x = "") +
  # stat_compare_means(comparisons = list(c("1:50", "1:100"), c("1:50", "1:200"), 
  #                                       c("1:50", "1:500"), c("1:100", "1:200"),
  #                                       c("1:100", "1:500"), c("1:200", "1:500"))) + 
  stat_compare_means(comparisons = list(c("1:50", "1:100"), c("1:50", "1:200"), 
                                        c("1:50", "1:500"))) + 
  theme(legend.position = "none")

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/SKBR3/MALBAC-scaling-correlation_bulk",
              height = 6, width = 7)

# Pairwise correlation
mb50_pw = cor(mb50_cn)
mb50_pw = data.table(sample = "1:50", 
                      V1 = rownames(mb50_pw)[row(mb50_pw)[upper.tri(mb50_pw, diag = F)]], 
                      V2 = colnames(mb50_pw)[col(mb50_pw)[upper.tri(mb50_pw, diag = F)]], 
                      pearson = c(mb50_pw[upper.tri(mb50_pw, diag = F)]))

mb100_pw = cor(mb100_cn)
mb100_pw = data.table(sample = "1:100", 
                      V1 = rownames(mb100_pw)[row(mb100_pw)[upper.tri(mb100_pw, diag = F)]], 
                      V2 = colnames(mb100_pw)[col(mb100_pw)[upper.tri(mb100_pw, diag = F)]], 
                      pearson = c(mb100_pw[upper.tri(mb100_pw, diag = F)]))

mb200_pw = cor(mb200_cn)
mb200_pw = data.table(sample = "1:200", 
                      V1 = rownames(mb200_pw)[row(mb200_pw)[upper.tri(mb200_pw, diag = F)]], 
                      V2 = colnames(mb200_pw)[col(mb200_pw)[upper.tri(mb200_pw, diag = F)]], 
                      pearson = c(mb200_pw[upper.tri(mb200_pw, diag = F)]))

mb500_pw = cor(mb500_cn)
mb500_pw = data.table(sample = "1:500", 
                     V1 = rownames(mb500_pw)[row(mb500_pw)[upper.tri(mb500_pw, diag = F)]], 
                     V2 = colnames(mb500_pw)[col(mb500_pw)[upper.tri(mb500_pw, diag = F)]], 
                     pearson = c(mb500_pw[upper.tri(mb500_pw, diag = F)]))

res_pw = rbindlist(list(mb50_pw, mb100_pw, mb200_pw, mb500_pw))

# Prepare for plotting
res_pw[, sample := factor(sample, levels = c("1:50", "1:100", "1:200", "1:500"))]
obs = res_pw[, .N, by = sample]

# Plotting
plt2 = ggplot(res_pw, aes(x = sample, y = pearson)) +
  geom_violin(aes(fill = sample, color = sample)) +
  geom_boxplot(width = .1) +
  geom_text(data = obs, aes(y = 0, label = paste("n =", N))) +
  scale_fill_viridis_d(begin = .4) +
  scale_color_viridis_d(begin = .4) +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, .2), labels = seq(0, 1, .2)) +
  labs(y = "Pairwise Pearson's Correlation", x = "") +
  stat_compare_means(comparisons = list(c("1:50", "1:100"), c("1:50", "1:200"), 
                                        c("1:50", "1:500"))) + 
  theme(legend.position = "none")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/SKBR3/MALBAC-scaling-pairwise-correlation",
              height = 6, width = 7)

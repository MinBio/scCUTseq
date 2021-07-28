## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script to calculate breadth of coverage

## Load/install packages
packages = c("data.table", "pbapply", "paletteer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

numthreads = 32

covhist = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/covhist/", full.names = T, recursive = T)

boc = pblapply(covhist, function(cov) {
  dt = fread(cov)

  # Calculate breadth of coverage
  boc = 1 - dt[V1 == "genome" & V2 == 0, V5]
  
  # Extract technique
  tech = tstrsplit(cov, "/")[[10]]
  
  # Extract sample
  sample = tstrsplit(cov, "/")[[11]]
  
  return(data.table(tech = tech,
                    sample = sample,
                    boc = boc))
}, cl = numthreads)

boc = rbindlist(boc)
boc[, tech := factor(tech, levels = c("scCUTseq", "DOP-PCR", "10X", "ACT"))]

plt = ggplot(boc, aes(x = tech, y = boc, group = sample)) +
  geom_quasirandom(size = .8, position = position_dodge(), dodge.width = .6, aes(color = tech)) +
  scale_color_paletteer_d("colorblindr::OkabeIto") +
  labs(x = "", y = "Breadth of Coverage") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/technical-validation/boc/breadth-of-coverage_techniques",
              height = 5, width = 7)

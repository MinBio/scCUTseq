## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "AneuFinder")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load models
models_NT = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO230+232+233/NZ84/MODELS/method-dnacopy/",
                    pattern = "5e\\+05", full.names = T)

models_APH = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO230+232+233/NZ85/MODELS/method-dnacopy/",
                        pattern = "5e\\+05", full.names = T)

# Select cells that meat quality thresholds
load(models_NT[1])
total_length = sum(width(model$bins))

nt_perc = unlist(lapply(models_NT, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) {
    (sum(width(model$segments[model$segments$copy.number != 2,])) / total_length) * 100    
  }
}))

aph_perc = unlist(lapply(models_APH, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) {
    (sum(width(model$segments[model$segments$copy.number != 2,])) / total_length) * 100
  }
}))

nt_seg = unlist(lapply(models_NT, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) {
    length(model$segments[model$segments$copy.number != 2,])
  }
}))

aph_seg = unlist(lapply(models_APH, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) {
    length(model$segments[model$segments$copy.number != 2,])
  }
}))

nt_cells = unlist(lapply(models_NT, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) {
    return(cell)
  }
}))

aph_cells = unlist(lapply(models_APH, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) {
    return(cell)
  }
}))


save_and_plot(heatmapGenomewide(aph_cells), "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO230+232+233/NZ85//PLOTS/method-dnacopy/hq-cells_5e5",
              width = 12, height = 5)
save_and_plot(heatmapGenomewide(nt_cells), "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO230+232+233/NZ84/PLOTS/method-dnacopy/hq-cells_5e5",
              width = 12, height = 5)

# plot %aneuploidy
aneu = data.table(variable = c(rep("72hr NT", length(nt_perc)), rep("72hr APH", length(aph_perc))),
                  aneuploidy = c(nt_perc, aph_perc))
aneu[, variable := factor(variable, levels = c("72hr NT", "72hr APH"))]

plt = ggplot(aneu, aes(x = variable, y = aneuploidy)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = variable)) +
  geom_jitter(width = 0.15) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Aneuploidy (%)",
       x = "") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/CNA-quantification/BICRO230+232+233_MCF10a-5e5-percentage-aneuploidy-APH72hr",
              width = 7, height = 7)

# Segments
aneu = data.table(variable = c(rep("72hr NT", length(nt_seg)), rep("72hr APH", length(aph_seg))),
                  aneuploidy = c(nt_seg, aph_seg))
aneu[, variable := factor(variable, levels = c("72hr NT", "72hr APH"))]

plt_seg = ggplot(aneu, aes(x = variable, y = aneuploidy)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = variable)) +
  geom_jitter(width = 0.15) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Aneuploidy (%)",
       x = "") +
  theme(legend.position = "none")

save_and_plot(plt_seg, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/CNA-quantification/BICRO230+232+233_MCF10a-5e5-segments-aneuploidy-APH72hr",
              width = 7, height = 7)

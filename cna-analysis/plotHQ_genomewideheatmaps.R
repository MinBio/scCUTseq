## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "AneuFinder")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/functions/plotGenomewideHeatmap.R")

models = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/BICRO237/NZ154/MODELS/method-dnacopy/",
                    pattern = "1e\\+06", full.names = T)

# Select cells that meat quality thresholds
models_hq = unlist(lapply(models, function(cell){
  load(cell)
  
  test_statements = c(model$qualityInfo$total.read.count > 3e5,
                      model$qualityInfo$avg.read.count > 15,
                      model$qualityInfo$spikiness < 0.7,
                      model$qualityInfo$num.segments < 100,
                      median(model$bins$copy.number) == 2)
  if(all(test_statements)) cell
}))


plt = plotGenomewideHeatmap(models_hq, 40)

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/NZ127-HQ-1mb",
              width = 18, height = 10)

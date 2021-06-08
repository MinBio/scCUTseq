## Author: Luuk Harbers
## Date: 2021-05-31
## Script for combining metrics and qualifiers

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

libs = c("MS81", "NZ169", "NZ170", "NZ186", "NZ250", "MS157")

cnvs = list.files("/mnt/AchTeraD/data/", full.names = T, recursive = T, pattern = "cnv.rds")
qual = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/cell_classifier/", full.names = T, pattern = ".csv")

cnvs = cnvs[grepl(paste(libs, collapse = "|"), cnvs)]
qual = qual[grepl(paste(libs, collapse = "|"), qual)]


metrics = lapply(cnvs, function(lib) {
  rds = readRDS(lib)
  metrics = rds$stats
  metrics[, library := gsub("/cnv.*|.*BICRO....", "", lib)]
  return(metrics)
})

metrics = rbindlist(metrics)

quals = lapply(qual, function(lib) {
  quality = fread(lib)
  quality[, library := gsub(".csv", "", basename(lib))]
  quality[, sample := gsub(".png", "", basename(sample))]
  return(quality)
})

quals = rbindlist(quals)

total = merge(metrics, quals, id.vars = c("sample", "library"))
write.table(total, "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/cell_classifier/metrics_LH.csv", sep = ",", col.names = T, row.names = F, quote = F)

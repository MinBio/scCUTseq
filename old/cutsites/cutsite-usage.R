## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "Rsamtools")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

cutsites = fread("/mnt/AchTeraD/Documents/Projects/COVseq/cutsite-distribution/sars-cov-2_NlaIII-MseI-cutsites.bed")


param = ScanBamParam(what = c("rname", "pos", "strand", "qwidth"))

files = list.files("/mnt/AchTeraD/data/BICRO248+249/TN11_S2/trimmed", pattern = ".bam$", full.names = T)
files = files[grepl("sample", files)]

res = lapply(files, function(bam) { 
  bam_data = as.data.frame(scanBam(bam, param = param))
  setDT(bam_data)
  setnames(bam_data, "rname", "chr")
  
  # setkeys
  setkey(bam_data, chr, pos)
  setkey(cutsites, V1, V2)
  
  # Check if cutsite has a read
  merged = cutsites[bam_data, roll = "nearest"]
  not_present = cutsites[!cutsites$V3 %in% merged$V3,]
  not_present[, sample := bam]
  return(not_present)
  })

total = rbindlist(res)
num = total[, .N, by = "sample"]
num[, performance := ifelse(grepl("sample1.bam|sample2.bam|sample3.bam|sample17.bam|sample21.bam|sample27.bam|sample28.bam|sample29.bam|sample30.bam", sample), "fail", "success")]

plt = ggplot(num, aes(x = performance, y = N, color = performance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  scale_color_viridis_d(end = 0.6)
  

occ = total[, .N, by = .(V1, V2, V3)]

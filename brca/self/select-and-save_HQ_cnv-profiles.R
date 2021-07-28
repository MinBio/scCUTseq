## Author: Luuk Harbers
## Date: 2021-12-07
## Script for selecting HQ profiles for phylogenetic analysis

## Load/install packages
packages = c("data.table", "pbapply")
sapply(packages, require, character.only = T)

# Set threads
nthreads = 32

# Select all libraries with breast cancer cells)
libraries = data.table(library = paste0("NZ", c(249:252, 255, 257)),
                       basepath = c(paste0("/mnt/AchTeraD/data/BICRO284/NZ", c(249:252, 255, 257), "/")))
# libraries = data.table(library = paste0("NZ", 169:174),
#                        basepath = c(paste0("/mnt/AchTeraD/data/BICRO277/NZ", 169:174, "/")))
# Loop through dirs and get copy numbers
total = pblapply(1:nrow(libraries), function(i) {
  rds = readRDS(paste0(libraries$basepath[i], "cnv/500000/cnv.rds"))
  bins = rds$bins
  stats = rds$stats
  dt = data.table(rds$copynumber)
  dt = dt[, colnames(dt) %in% stats[stats$classifier_prediction == "good", sample], with = F]
  
  # Set colnames to include library
  setnames(dt, paste0(libraries$library[i], "_", colnames(dt)))
  
  # Return dt
  return(cbind(bins, dt))
}, cl = nthreads)

merged = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), total)
merged = merged[complete.cases(merged)]

write.table(merged[, 4:ncol(merged)], "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca_2-integer_cn-HQ.tsv",
            quote = F, row.names = F, col.names = T, sep = "\t")



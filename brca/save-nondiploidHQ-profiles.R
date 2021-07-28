## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for saving rds of non-diploid HQ profiles

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
  # dt = dt[, colMeans(dt) > 2.15 & colMeans(dt) < 2.5, with = F] #BRCA1
  #dt = dt[, colMeans(dt) > 2.05 | colMeans(dt) < 1.95, with = F] #BRCA2
  # Set colnames to include library
  setnames(dt, paste0(libraries$library[i], "_", colnames(dt)))
  
  # Return dt
  return(cbind(bins, dt))
}, cl = nthreads)

merged = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), total)
merged = merged[complete.cases(merged)]

# # Extra selection ONLY BRCA 2
# hc = hclust(dist(t(merged[, 4:ncol(merged)])), method = "average")
# select = cutree(hc, h = 46)
# merged = merged[, c(T, T, T, colnames(merged[, 4:ncol(merged)]) %in% names(select[select == 1])), with = F]


saveRDS(merged, "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca2_nondiploid-HQ.rds")

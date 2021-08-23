## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)

nthreads = 32

annot = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/P6.tsv", header = F)
# Get file location
files = list.files("/mnt/AchTeraD/data", recursive = T, full.names = T, pattern = "cnv.rds")
files = files[grepl(paste(annot$V1, collapse = "|"), files)]
files_dt = data.table(files = files,
                      library = gsub("\\/cnv.*|.*BICRO....", "", files))

annot = merge(annot, files_dt, by.x = "V1", by.y = "library")
setnames(annot, c("library", "section", "file"))

# Loop through dirs and get copy numbers
total = pblapply(1:nrow(annot), function(i) {
  rds = readRDS(annot[i, file])
  bins = rds$bins
  stats = rds$stats
  dt = data.table(rds$copynumber)
  dt = dt[, colnames(dt) %in% stats[stats$classifier_prediction == "good", sample], with = F]
 
  # Set colnames to include library
  setnames(dt, paste0(annot$library[i], "_", colnames(dt)))
  
  # Return dt
  return(cbind(bins, dt))
}, cl = nthreads)

merged = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), total)
merged = merged[complete.cases(merged)]

saveRDS(merged, "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/prostate_p6.rds")

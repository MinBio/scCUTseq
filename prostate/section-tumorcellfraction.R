## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting TCF

## Load/install packages
packages = c("data.table", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")
source("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/plotHeatmap.R")

annot = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/P3.tsv", header = F)
# Get file location
files = list.files("/mnt/AchTeraD/data", recursive = T, full.names = T, pattern = "cnv.rds")
files = files[grepl(paste(annot$V1, collapse = "|"), files)]
files_dt = data.table(files = files,
                      library = gsub("\\/cnv.*|.*BICRO....", "", files))

annot = merge(annot, files_dt, by.x = "V1", by.y = "library")
setnames(annot, c("library", "section", "file"))

# Read in copynumber data and get fraction of tumor cells vs diploid cells
res = pblapply(1:nrow(annot), function(i) {
  rds = readRDS(annot[i, file])
  bins = rds$bins
  stats = rds$stats
  dt = data.table(rds$copynumber)
  dt = dt[, colnames(dt) %in% stats[stats$classifier_prediction == "good", sample], with = F]
  #dt = dt[, colMeans(dt) < 1.94 | colMeans(dt) > 1.96, with = F]
  
  # Keep nondiploid profiles
  nondiploid = apply(dt[which(bins$chr != "X")], 2, function(x) sum(x != 2))
  alt_x = apply(dt[which(bins$chr == "X")], 2, function(x) sum(x != 1))
  
  # Keep
  keep = unique(c(names(nondiploid[nondiploid > 0]), names(alt_x[alt_x > 0])))
  
  # Return dt
  return(data.table(section = annot[i, section],
                    nondiploid = length(keep),
                    n = ncol(dt),
                    nondiploid_perc = length(keep)/ncol(dt)))
})

total = rbindlist(res)

# Get coordinates
total[, x := as.numeric(gsub("L|C.", "", section))]
total[, y := as.numeric(gsub("L.|C", "", section))]

# Fill NAs
total = complete(total, expand(total, x, y))

plt = ggplot(total, aes(x = x, y = y, fill = nondiploid_perc)) +
  geom_tile() +
  scale_fill_distiller(name = "Tumor Cell\nFraction", palette = "Reds", direction = 1, na.value = "grey") +
  geom_hline(yintercept = seq(from = .5, to = max(total$y), by = 1)) +
  geom_vline(xintercept = seq(from = .5, to = max(total$x), by = 1)) +
  scale_y_reverse(expand = c(0, 0 ), breaks = seq(1, max(total$y)), labels = seq(1, max(total$y))) + 
  scale_x_reverse(expand = c(0, 0 ), breaks = seq(1, max(total$x)), labels = seq(1, max(total$x))) +
  theme(axis.title = element_blank())

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/TCF_blocks/P3_TCF",
              height = 4, width = 8)  

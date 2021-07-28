packages = c("data.table", "ggplot2", "umap", "Rtsne")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/"

rda = list.files(base_path, pattern = "-250000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-250000.Rds", recursive = T, full.names = T)

rda = rda[grepl("NZ39|NZ40|BICRO243|BICRO230\\+232\\+233/NZ84", rda)]
rds = rds[grepl("NZ39|NZ40|BICRO243|BICRO230\\+232\\+233/NZ84", rds)]

#types = c("SKBR3_live", "SKBR3_fixed", "MCF10A_fixed", "SKBR3-in-IMR90_fixed", rep("TK6_fixed", 16))
types = c("SKBR3_live", "SKBR3_fixed", "MCF10A_fixed", "SKBR3-in-IMR90_fixed")

data = lapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  return(list(data.table(psoptim), stats))
})

# Set thresholds
min_reads = 2e6
max_spikiness = 0.5
min_avgreads = 50

total = lapply(1:length(data), function(lib) {
  # Get usable cells and filter dt
  usable_cells = data[[lib]][[2]][reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = data[[lib]][[1]][, usable_cells, with = F]
  setnames(dt, paste0(types[lib], "-", colnames(dt)))
  
  if(lib > 2){
    # Select diploid cells
    dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
  }
  
  # return dt
  return(dt)
})

total = do.call(cbind, total)
total = total[1:nrow(total)-1,]

# Get annotation
annot = data.table(sample = colnames(total), group = gsub("SKBR3-in-|-.*", "", colnames(total)))

# Get defaults to change
new_config = umap.defaults
new_config$n_neighbors = 15


total_umap = umap(t(total), config = new_config, method = "umap-learn")
umap_dt = data.table(x = total_umap$layout[, 1],
                    y = total_umap$layout[, 2],
                    group = annot$group)

plot = ggplot(umap_dt, aes(x = x, y = y, color = group)) +
  geom_point(size = 2) +
  scale_color_viridis_d(begin = 0.1, direction = 1) +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Cell type")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/celltypes_heatmap/umap_clustering-celltypes_250kb",
              height = 7, width = 7)


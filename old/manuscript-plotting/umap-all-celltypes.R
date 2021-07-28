packages = c("data.table", "ggplot2", "umap", "Rtsne")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/"

rda = list.files(base_path, pattern = "-250000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-250000.Rds", recursive = T, full.names = T)

rda = rda[grepl("NZ39|NZ40|BICRO243|BICRO230\\+232\\+233/NZ84|P18158_combined", rda)]
rds = rds[grepl("NZ39|NZ40|BICRO243|BICRO230\\+232\\+233/NZ84|P18158_combined", rds)]

types = c("SKBR3", "SKBR3", "MCF10A", "IMR90", rep("TK6", 16))

data = lapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  return(list(data.table(psoptim), stats))
})

# Set thresholds
min_reads = 1e3
max_spikiness = 0.55
min_avgreads = 50

total = lapply(1:length(data), function(lib) {
  # Get usable cells and filter dt
  usable_cells = data[[lib]][[2]][reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = data[[lib]][[1]][, usable_cells, with = F]
  setnames(dt, paste0(types[lib], "-", colnames(dt)))
  
  if(lib > 4) {
    startdel=6770
    enddel=6797
    
    outsidestart=6750
    outsideend=6817
    
    lapply(colnames(dt), function(cell) {
      inside = sum(dt[, ..cell][startdel:enddel] < 2)
      outside = sum(dt[, ..cell][c(outsidestart:startdel, enddel:outsideend)] < 2)
      if(inside > 0.6 * (enddel - startdel) & outside < 0.3 * (startdel - outsidestart + outsideend - enddel)){
        setnames(dt, cell, paste0("CAS9 deleted ", cell))
      }
    })
  }
  
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
new_config$n_neighbors = 38
new_config$min_dist = 0.8


total_umap = umap(t(total), config = new_config, method = "umap-learn")
umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     group = annot$group)
umap_dt[, group := factor(group, levels = c("TK6", "CAS9 deleted TK6", "IMR90", "MCF10A", "SKBR3"))]

plot = ggplot(umap_dt, aes(x = x, y = y, color = group)) +
  geom_point(size = 1.5) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Cell type")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/celltypes_heatmap/umap_clustering-all-celltypes_250kb",
              height = 7, width = 7)


## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "ggplot2", "umap", "pbapply", "ggdendro", "ggalt")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 16
set.seed(666)

base_path = "/mnt/AchTeraD/data/"

rda = list.files(base_path, pattern = "-500000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-500000.Rds", recursive = T, full.names = T)
rda = rda[grepl("BICRO241|BICRO249", rda)]
rds = rds[grepl("BICRO241|BICRO249", rds)]
annot = list(b241 = fread("/mnt/AchTeraD/data/BICRO241/BICRO241_384 Barcodes.tsv", header = T),
             b249 = fread("/mnt/AchTeraD/data/BICRO249/BICRO249_384 Barcodes.tsv", header = T))

runs = c("BICRO241", "BICRO249")

# Set QC thresholds
min_reads = 3e5
max_spikiness = 0.55
min_avgreads = 50

total = pblapply(1:length(rda), function(lib) {
  load(rda[lib])
  info = readRDS(rds[lib])
  stats = info$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  # Filter out cells from other sequencing
  if(lib == 1){
    brca_cells = annot[[lib]][cell_type == "Bca biopsy", BC_original]
    dt = dt[, colnames(dt) %in% brca_cells, with = F]
  }

  # Filter out low quality cells
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, colnames(dt) %in% usable_cells, with = F]
  
  # Filter out low quality cells that still manage to sneak through previous filtering
  colmeans = colMeans(dt)
  dt = dt[, colmeans < 5, with = F]
  
  # Set name with run info
  setnames(dt, paste0(names(annot)[lib], "-", colnames(dt)))
  
  # Merge with bins
  dt = cbind(info$binbed[[1]], dt)
  
  # Return dt
  return(dt)
}, cl = nthreads)

# Bind total and get annotations
total = merge(total[[1]], total[[2]], by = c("chr", "start" , "end"), sort = F)

# Make umap
config = umap.defaults
config$n_neighbors = 12
config$min_dist = 0.2

total_umap = umap(t(total[, 4:ncol(total)]), method = "umap-learn", config = config)

umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     sample = colnames(total[, 4:ncol(total)]))

plot = ggplot(umap_dt, aes(x = x, y = y)) +
  geom_point(size = 2) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2")



# DEFINE CLUSTERS -- MANUALLY
umap_dt[x < 0, clusters := "Cluster-1"]
umap_dt[x > 0 & x < 4 & y > 3, clusters := "Cluster-2"]
umap_dt[x < 10 & x > 4 & y > 2.6 & y < 5, clusters := "Cluster-3"]
umap_dt[x > 0 & x < 10 & y > -3 & y < 2.6, clusters := "Cluster-4"]
umap_dt[x > 15, clusters := "Cluster-5"]

# Plot with circles
plot = ggplot(umap_dt, aes(x = x, y = y, color = clusters)) +
  geom_point(size = 2) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Block")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/bc241_249-UMAP-clusters",
              height = 6, width = 9)

# Prepare for plotting heatmap
bins = total[, 1:3]

# Get cumulative locations
bins[, bin := seq_along(chr)]
bins[, end_cum := cumsum((end - start) + 1)]
bins[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Make chr_bounds
chr_bounds = bins[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
chr_bounds = chr_bounds %>% 
  mutate(mid = round(min + (max-min) / 2,0),
         end_bp=cumsum(as.numeric(chrlen_bp)), 
         start_bp = end_bp - chrlen_bp, 
         mid_bp = round((chrlen_bp / 2) + start_bp, 0))

#Colors
colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# sample 8 profiles from each cluster
#umap_dt = umap_dt[, .SD[sample(.N, min(.N, 4))], by = clusters]
cols = umap_dt$sample

# Set factors
setorder(umap_dt, clusters)
umap_dt[, sample := factor(sample, levels = sample)]

# Plot annotation bar
annot_plt = ggplot(umap_dt, aes(x = 1, y = sample, fill = clusters)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme_void() +
  theme(legend.position = "left")

dt = cbind(bins, total[, ..cols])

dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

# Set sample order
dt_melt[, variable := factor(variable, levels = umap_dt$sample)]

# Plot heatmap
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 2) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

combined = cowplot::plot_grid(annot_plt,
                              heatmap,
                              align = "h", rel_widths = c(0.08, 1), ncol = 2)

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/bc241_249-heatmap-clusters",
              height = 12, width = 20)

## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting HQ heatmap if selected libraries

## Load/install packages
packages = c("data.table", "ggplot2", "ggdendro", "matrixStats")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/"

rda = list.files(base_path, pattern = "-250000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-250000.Rds", recursive = T, full.names = T)

rda = rda[grepl("NZ39|NZ40|BICRO243|BICRO230\\+232\\+233/NZ84", rda)]
rds = rds[grepl("NZ39|NZ40|BICRO243|BICRO230\\+232\\+233/NZ84", rds)]

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
max_spikiness = 0.50
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

# Prepare for plotting heatmap
stats = readRDS(rds[1])
bins = stats$binbed[[1]]
setDT(bins)

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
colors = c("#496bab", "#9fbdd7", "#c1c1c1", "#e9c47e", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

total = do.call(cbind, total)
dt = data.table(cbind(bins, total))
dt = dt[1:nrow(dt)-1,]

# Clustering
hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
dhc = as.dendrogram(hc)

# Rectangular lines
ddata = dendro_data(dhc, type = "rectangle")

# Plot Dendrogram
dendro = ggplot(ggdendro::segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.004, 0.004)) +
  theme_dendro()

# Prepare for heatmap
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

# # Set sample order
dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

# Plot annotation
annot = data.table(samples = factor(colnames(dt[, 7:ncol(dt)]), levels = ddata$labels$label),
                   cell_type = gsub("SKBR3-in-|-.*", "", colnames(dt[, 7:ncol(dt)])))
annot[, cell_type := factor(cell_type, levels = c("IMR90_fixed", "MCF10A_fixed", "SKBR3_fixed", "SKBR3_live"))]


annot_plt = ggplot(annot, aes(x = 1, y = samples, fill = cell_type)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d(begin = 0.4) +
  theme_void() 

# Plot heatmap
heatmap = ggplot(dt_melt) +
  ggrastr::rasterise(geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 2),
                     dpi = 2400) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

# Combine plots and save
combined_noleg = cowplot::plot_grid(dendro,
                                    annot_plt + theme(legend.position = ""),
                                    heatmap,
                                    align = "h", rel_widths = c(0.1, 0.02, 1), ncol = 3)

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = cowplot::plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.05))

cairo_ps("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/celltypes_heatmap/2mil+250kb-heatmap_celltypes_cairo_ps.eps", 
         onefile = TRUE, height=6, width=16, family="Helvetica", pointsize=8, antialias="none", fallback_resolution = 2400)
print(combined)
dev.off()

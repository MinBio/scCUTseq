## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting HQ heatmap if selected libraries

## Load/install packages
packages = c("data.table", "ggplot2", "ggdendro")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/BICRO246"

rda = list.files(base_path, pattern = "-250000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-250000.Rds", recursive = T, full.names = T)
annot = fread("/mnt/AchTeraD/data/BICRO246/annot-barcodes.tsv", header = F)

load(rda)
stats = readRDS(rds)
stats = stats$stats[[1]]
setDT(stats)
cn = data.table(psoptim)

# Set thresholds
min_reads = 3e5
max_spikiness = 0.45
min_avgreads = 50

# Get usable cells and filter dt
usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
dt = cn[, usable_cells, with = F]

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

dt = data.table(cbind(bins, dt))

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
used = data.table(V1 = colnames(dt[, 7:ncol(dt)]))
annot = merge(used, annot, by = "V1", all.x = TRUE)
annot[, V3 := factor(V3, levels = c("Patient2", "Patient3", "Patient4"))]
annot = annot[match(ddata$labels$label, annot$V1),]
annot[, V1 := factor(V1, levels = V1)]

annot_plt = ggplot(annot, aes(x = 1, y = V1, fill = V3)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d(begin = 0.4) +
  theme_void() 

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

# Combine plots and save
combined_noleg = cowplot::plot_grid(dendro,
                                    annot_plt + theme(legend.position = ""),
                                    heatmap,
                                    align = "h", rel_widths = c(0.1, 0.02, 1), ncol = 3)

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = cowplot::plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.05))

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/genomewideheatmap-BICRO246_HQ-annotated-500kb",
              height = 4, width = 16)

config = umap.defaults
config$n_neighbors = 5


total_umap = umap(t(dt[, 7:ncol(dt)]), method = "umap-learn", config = config)

umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     group = annot$V3,
                     sample = annot$V1)

plot = ggplot(umap_dt, aes(x = x, y = y, color = group)) +
  geom_point(size = 2) +
  scale_color_viridis_d(begin = 0.4) +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Patient")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/UMAP-BICRO246_HQ-annotated-500kb",
              height = 6, width = 7.5)

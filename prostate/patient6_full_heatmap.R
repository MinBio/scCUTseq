## Author: Luuk Harbers
## Date: 2021-05-14
## Script for heatmap generation of patient 3

## Load/install packages
packages = c("data.table", "ggplot2", "naturalsort", "pbapply", "ggdendro", "ggalt", "dendextend", "umap")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 32

# Block annotation
blocks = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/patient6_blocks.tsv")

base_path = "/mnt/AchTeraD/data/BICRO284/"

rda = list.files(base_path, pattern = "-500000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-500000.Rds", recursive = T, full.names = T)
rda = rda[grepl(paste(blocks$library, collapse = "|"), rda)]
rds = rds[grepl(paste(blocks$library, collapse = "|"), rds)]
# Set QC thresholds
min_reads = 3e5
max_spikiness = 0.5
min_avgreads = 50

total = pblapply(1:length(rda), function(lib) {
  load(rda[lib])
  rds = readRDS(rds[lib])
  bins = rds$binbed[[1]]
  setDT(bins)
  stats = rds$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  # Filter out low quality cells
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, colnames(dt) %in% usable_cells, with = F]
  
  # # Filter out low quality cells that still manage to sneak through previous filtering
  # colmeans = colMeans(dt)
  # dt = dt[, colmeans < 5, with = F]
  
  # Set name with block info
  currentlib = gsub("\\/out.*", "", rda[lib])
  currentlib = gsub(".*\\/", "", currentlib)
  
  setnames(dt, paste0(blocks[library == currentlib]$block, "_", colnames(dt)))
  
  # Return dt
  return(cbind(bins, dt))
}, cl = nthreads)

# Bind total and get annotations
total = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end")), total)
total = na.omit(total)

# Prepare for plotting heatmap
total = total[naturalorder(chr)]

# Get cumulative locations
total[, bin := seq_along(chr)]
total[, end_cum := cumsum((end - start) + 1)]
total[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Make chr_bounds
chr_bounds = total[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
chr_bounds = chr_bounds %>% 
  mutate(mid = round(min + (max-min) / 2,0),
         end_bp=cumsum(as.numeric(chrlen_bp)), 
         start_bp = end_bp - chrlen_bp, 
         mid_bp = round((chrlen_bp / 2) + start_bp, 0))

#Colors
colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# Make dendrogram
hc = hclust(dist(t(total[, .SD, .SDcols = !c("chr", "start", "end", "bin", "end_cum", "start_cum")])), method = "average")
# dhc = as.dendrogram(hc)
# 
# # Rectangular lines
# ddata = dendro_data(dhc, type = "rectangle")
# 
# # Plot Dendrogram
# dendro = ggplot(ggdendro::segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
#   coord_flip() + 
#   scale_y_reverse(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0.0, 0.0)) +
#   theme_dendro()

# Set order
sample_order = colnames(total[, .SD, .SDcols = !c("chr", "start", "end", "bin", "end_cum", "start_cum")])[order.hclust(hc)]

# Plot annotation bar
# blocks[, sample := factor(sample, levels = ddata$labels$label)]
blocks = data.table(sample = sample_order, block = gsub("_.*", "", sample_order))
blocks[, sample := factor(sample, levels = sample_order)]

# Make umap 
config = umap.defaults
config$n_neighbors = 25 
config$min_dist = 0.2

total_umap = umap(t(total[, .SD, .SDcols = !c("chr", "start", "end", "bin", "end_cum", "start_cum")]), 
                  config = config)

umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     block = blocks$block,
                     sample = blocks$sample)

plot = ggplot(umap_dt, aes(x = x, y = y, color = block)) +
  geom_point(size = 1.5) +
  scale_color_hue() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Block", title = paste0("n = ", length(sample_order)))

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/UMAP_patient6",
              width = 9, height = 7)

# Plot genome profile heatmap
annot_plt = ggplot(blocks, aes(x = 1, y = sample, fill = block)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_hue() +
  theme_void() +
  theme(legend.position = "none")

dt_melt = melt(total, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
# dt_melt[, variable := factor(variable, levels = ddata$labels$label)]
dt_melt[, variable := factor(variable, levels = sample_order)]

# Plot heatmap
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = .2) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number", title = paste0("n = ", length(sample_order))) + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

# Combine plots and save
# combined_noleg = cowplot::plot_grid(dendro,
#                                     annot_plt + theme(legend.position = ""),
#                                     heatmap,
#                                     align = "h", axis = "tb", rel_widths = c(0.15, 0.01, 1), ncol = 3)
# 
# legend = cowplot::get_legend(
#   annot_plt + guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom"))
# 
# combined = cowplot::plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.05))
combined_noleg = cowplot::plot_grid(annot_plt + theme(legend.position = ""),
                                    heatmap,
                                    align = "h", axis = "tb", rel_widths = c(0.02, 1), ncol = 2)

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = cowplot::plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.1))
save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/patient6_genomewideheatmap-novaseq",
              height = 20, width = 20)

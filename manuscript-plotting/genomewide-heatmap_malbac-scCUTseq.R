## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting heatmap of selected cells

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/"
runs = c("BICRO218", "BICRO221")

files = unlist(sapply(runs, function(run) list.files(paste0(base_path, "/", run), recursive = T, pattern = ".Rd*", full.names = T)))
files = files[grepl("BICRO221|NZ39|NZ40", files)]

rda = files[grepl("Rda", files)]
rds = files[grepl("Rds", files)]

# Set thresholds
min_reads = 1e6
max_spikiness = 0.5
min_avgreads = 100

# Loop through file and select usable cells
data = lapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]

  # Select cells
  setDT(stats)
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  
  psoptim = data.table(psoptim)
  psoptim = psoptim[, usable_cells, with = F]
  setnames(psoptim, paste(names(rds[lib]), colnames(psoptim), sep = "_"))
  
  return(psoptim)
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

# Combine data.tables
total = bind_cols(data)

# combine data
dt = data.table(cbind(bins, total))

# Make dendrogram
hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
dhc = as.dendrogram(hc)

# Rectangular lines
ddata = dendro_data(dhc, type = "rectangle")

# Plot Dendrogram
dendro = ggplot(ggdendro::segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() + 
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.018, 0.018)) +
  theme_dendro()

# Plot annotation bar
annot = data.table(cells = colnames(total), 
                   type = c(rep("scCUTseq - live", 6), 
                            rep("scCUTseq - fixed", 14),
                            "MALBAC_1:200 - fixed", "MALBAC_1:200 - live",
                            "MALBAC_1:200 - fixed", "MALBAC_1:200 - live",
                            "MALBAC_1:200 - live", "MALBAC_1:200 - fixed"))
annot[, cells := factor(cells, levels = ddata$labels$label)]
annot_plt = ggplot(annot, aes(x = 1, y = cells, fill = type)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d() +
  theme_void() 

# Prepare for heatmap
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

# Set sample order
dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

# Plot heatmap 20K+
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 5) +
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
                                    align = "h", rel_widths = c(0.15, 0.01, 1), ncol = 3)

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = cowplot::plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.05))

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/malbac-scCUTseq/genomewideheatmap-malbac_scCUTseq",
              height = 5, width = 16)

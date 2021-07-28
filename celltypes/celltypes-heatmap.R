## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")
source("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/plotHeatmap.R")

# Read in cnv.rds and extract HQ profiles
rds_files = data.table(celltype = c("SKBR3_live", "SKBR3_fixed", "IMR90", "MCF10A"),
                       file = c("/mnt/AchTeraD/data/BICRO218/NZ39/cnv/500000/cnv.rds", 
                                "/mnt/AchTeraD/data/BICRO218/NZ40/cnv/500000/cnv.rds",
                                "/mnt/AchTeraD/data/BICRO231/NZ86/cnv/500000/cnv.rds", 
                                "/mnt/AchTeraD/data/BICRO232/NZ84/cnv/500000/cnv.rds"))

cnv = pblapply(1:nrow(rds_files), function(i) {
  rds = readRDS(rds_files$file[i])
  cn = rds$copynumber[, rds$stats[classifier_prediction == "good", sample], with = F]
  setnames(cn, paste0(rds_files$celltype[i], "-", colnames(cn)))
  cn = cbind(rds$bins, cn)
})

# Bind cols
cn = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), cnv)
cn = cn[complete.cases(cn)]

# heatmap = plotHeatmap(cn[, 4:ncol(cn)], cn[, 1:3], linesize = 3, rasterize = T)
## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")
source("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/plotHeatmap.R")

# Read in cnv.rds and extract HQ profiles
rds_files = data.table(celltype = c("SKBR3_live", "SKBR3_fixed", "IMR90", "MCF10A"),
                       file = c("/mnt/AchTeraD/data/BICRO218/NZ39/cnv/500000/cnv.rds", 
                                "/mnt/AchTeraD/data/BICRO218/NZ40/cnv/500000/cnv.rds",
                                "/mnt/AchTeraD/data/BICRO231/NZ86/cnv/500000/cnv.rds", 
                                "/mnt/AchTeraD/data/BICRO232/NZ84/cnv/500000/cnv.rds"))

cnv = pblapply(1:nrow(rds_files), function(i) {
  rds = readRDS(rds_files$file[i])
  cn = rds$copynumber[, rds$stats[classifier_prediction == "good", sample], with = F]
  setnames(cn, paste0(rds_files$celltype[i], "-", colnames(cn)))
  cn = cbind(rds$bins, cn)
})

# Bind cols
cn = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), cnv)
cn = cn[complete.cases(cn)]

# heatmap = plotHeatmap(cn[, 4:ncol(cn)], cn[, 1:3], linesize = 3, rasterize = T)

# Get cumulative locations
cn[, bin := seq_along(chr)]
cn[, end_cum := cumsum((end - start) + 1)]
cn[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Make chr_bounds
chr_bounds = cn[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
chr_bounds = chr_bounds %>% 
  mutate(mid = round(min + (max-min) / 2,0),
         end_bp=cumsum(as.numeric(chrlen_bp)), 
         start_bp = end_bp - chrlen_bp, 
         mid_bp = round((chrlen_bp / 2) + start_bp, 0))

#Colors
colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# Clustering
dt = cn[, c(1:3, (ncol(cn)-2):ncol(cn), 4:(ncol(cn)-3)), with = F]
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
                   cell_type = gsub("-.*", "", colnames(dt[, 7:ncol(dt)])))

annot_plt = ggplot(annot, aes(x = 1, y = samples, fill = cell_type)) +
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

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/cell-types/imr90-mcf10a-skbr3-heatmap",
              width = 28, height = 9)

## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "pbapply", "ggdendro")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 32

lib_dirs_1 = list.dirs("/mnt/AchTeraD/data/BICRO277", recursive = F)
lib_dirs_1 = lib_dirs[grepl(paste(c("MS77", "MS78", "MS79", "MS80", "MS81", "NZ185", "NZ186", "NZ187", "NZ188", "NZ189"), collapse = "|"), lib_dirs_1)]
lib_dirs_2 = list.dirs("/mnt/AchTeraD/data/BICRO284", recursive = F)
lib_dirs_2 = lib_dirs[grepl(paste(c("MS101", "MS143", "MS144", "NZ208", "NZ212", "NZ215", "NZ216"), collapse = "|"), lib_dirs_2)]

lib_dirs = c(lib_dirs_1, lib_dirs_2, list.dirs("/mnt/AchTeraD/data/BICRO278", recursive = F))


cn = pblapply(lib_dirs, function(lib) {
  rds = readRDS(paste0(lib, "/cnv/500000/cnv.rds"))
  cn = rds$copynumber[, rds$stats[classifier_prediction == "good", sample], with = F]
  setnames(cn, paste0(basename(lib), "_", colnames(cn)))
  cn = cbind(rds$bins, cn)
}, cl = nthreads)

# Bind cols
merged = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), cn)
merged = merged[complete.cases(merged)]

# Select bins
bins = merged[, 1:3]

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

# Make dt
dt = cbind(bins, merged[, 4:ncol(merged)])

# Distance and clustering
dist = dist(t(dt[, 7:ncol(dt)]))
hc = hclust(dist, method = "average")

# Make dendrogram
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

# Set sample and chromosome orderorder
dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

# Plot heatmap
heatmap = ggplot(dt_melt) +
  ggrastr::rasterize(geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = .5)) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number", subtitle = paste0("n = ", ncol(dt) -7)) + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

plt = plot_grid(dendro, heatmap,  align = "h", rel_widths = c(.3, 2), ncol = 2)
save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/BRCA1/genomewideheatmap-HQ-250kb",
              width = 24, height = 12)

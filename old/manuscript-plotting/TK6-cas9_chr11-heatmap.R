## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting TK6 cells deletion

## Load/install packages
packages = c("data.table", "ggplot2", "ggdendro")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

base_path = "/mnt/AchTeraD/data/ngi/P18158_combined/"

rda = list.files(base_path, pattern = "-250000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-250000.Rds", recursive = T, full.names = T)
libraries = basename(list.dirs(base_path, recursive = F))
libraries = libraries[!grepl("fastq|fixed", libraries)]

data = lapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  return(list(data.table(psoptim), stats))
})

# Set thresholds
min_reads = 3e5
max_spikiness = 0.55
min_avgreads = 50
diploid = TRUE


total = lapply(1:length(data), function(lib) {
  # Get usable cells and filter dt
  usable_cells = data[[lib]][[2]][reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = data[[lib]][[1]][, usable_cells, with = F]
  
  # Select majorly diploid
  if(diploid) dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
  
  setnames(dt, paste0(libraries[lib], "_", colnames(dt)))
  
  # return dt
  return(dt)
})
total = do.call(cbind, total)


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

dt = data.table(cbind(bins, total))
# Select diploid cells
dt = dt[, c(rep(TRUE, 6), apply(dt[, 7:ncol(dt)], 2, function(x) mean(x) > 1.9 & mean(x) < 2.1)), with = F]

# Set start and end bin for chr11:118,359,229-125,769,283
delstart = 6770
delend = 6797
plotstart = 6743
plotend = 6824

dels = apply(dt[, 7:ncol(dt)], 2, function(x) {
  sum(x[delstart:delend] < 2)
})

dt = dt[plotstart:plotend, c("chr", "start", "end", "bin", "end_cum", "start_cum", names(dels[dels > 13])), with = F]

# Order
hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
dhc = as.dendrogram(hc)

# Rectangular lines
ddata = dendro_data(dhc, type = "rectangle")

# Prepare for heatmap
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

dt_melt[, variable := factor(variable, levels = ddata$labels$label)]
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 2) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), breaks = c(1771143146, 1778159884, 1785485500, 1792486784),
                     labels = c("-7mb", "gRNA 1", "gRNA2", "+7mb")) +
  geom_hline(yintercept = c(1778159884, 1785485500), linetype = 2, color = "red") +
  theme(axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

save_and_plot(heatmap, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-heatmaps/chr11-deletion_zoomin",
              height = 4, width = 6)


# cairo_ps("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-heatmaps/chr11-deletion_zoomin_rast_cairo_ps.eps",
#          onefile = TRUE, height=6.2, width=8, family="Helvetica", pointsize=8, antialias="none",
#          fallback_resolution = 2400)
# heatmap
# dev.off()

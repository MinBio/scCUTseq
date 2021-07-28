## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for TK6 analysis

## Load/install packages
packages = c("data.table", "pbapply", "uwot", "ggdendro", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
numthreads = 32
set.seed(1337)
binsize = "250000"

treated_dirs = list.dirs("/mnt/AchTeraD/data/BICRO276", recursive = F)
untreated_dirs = list.dirs("/mnt/AchTeraD/data/BICRO277", recursive = F)
untreated_dirs = untreated_dirs[grepl("NZ175|NZ179|NZ180|NZ181|NZ182|NZ183", untreated_dirs)]
bins = readRDS(paste0(treated_dirs[1], "/cnv/", binsize, "/cnv.rds"))$bins
setkey(bins, chr, start, end)

# Load in cnv info
treated = pblapply(treated_dirs, function(lib) {
  cnv = readRDS(paste0(lib, "/cnv/", binsize, "/cnv.rds"))
  
  hq_cells = cnv$stats[classifier_prediction == "good", sample]
  cn = cbind(cnv$bins, cnv$copynumber[, ..hq_cells])
  setnames(cn, c("chr", "start", "end", paste0(basename(lib), "_", hq_cells)))
  return(cn)
})

treated = Reduce(function(x, y) merge(x = x, y = y, by = c("chr", "start", "end")), treated)

untreated = pblapply(untreated_dirs, function(lib) {
  cnv = readRDS(paste0(lib, "/cnv/", binsize, "/cnv.rds"))
  
  hq_cells = cnv$stats[classifier_prediction == "good", sample]
  cn = cbind(cnv$bins, cnv$copynumber[, ..hq_cells])
  setnames(cn, c("chr", "start", "end", paste0(basename(lib), "_", hq_cells)))
  return(cn)
})

untreated = Reduce(function(x, y) merge(x = x, y = y, by = c("chr", "start", "end")), untreated)

# Make umap of all cells
dt = merge(treated, untreated)

# # Get defaults to change
# new_config = umap.defaults
# new_config$n_neighbors = 25
# 
# total_umap = umap(t(dt[, 4:ncol(dt)]), config = new_config)
# umap_dt = data.table(x = total_umap$layout[, 1],
#                      y = total_umap$layout[, 2],
#                      sample = rownames(total_umap$layout),
#                      group = ifelse(rownames(total_umap$layout) %in% colnames(treated), "Treated", "Untreated"))
# 
# plot = ggplot(umap_dt, aes(x = x, y = y, color = group)) +
#   geom_point(size = 2) +
#   scale_color_npg() +
#   labs(x  = "UMAP 1", y = "UMAP 2", color = "")
# 
# save_and_plot(plot, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/umap-", binsize),
#               width = 7, height = 7)

# UMAP of chr 11
# Get defaults to change
# new_config = umap.defaults
# new_config$n_neighbors = 25
# 
# total_umap = umap(t(dt[chr == "11", 4:ncol(dt)]), config = new_config)
set.seed(45)
total_umap = umap(t(dt[chr == "11", 4:ncol(dt)]), n_neighbors = 24, spread = 1, min_dist = 0)
umap_dt = data.table(x = total_umap[, 1],
                     y = total_umap[, 2],
                     sample = colnames(dt[, 4:ncol(dt)]),
                     group = ifelse(colnames(dt[, 4:ncol(dt)]) %in% colnames(treated), "Treated", "Untreated"))

plot = ggplot(umap_dt, aes(x = x, y = y, color = group)) +
  geom_point(size = 1.5) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "")

save_and_plot(plot, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/chr11_umap-", binsize),
              width = 7, height = 7)


# Set clones
clones = dbscan::hdbscan(umap_dt[,c(1:2)],
                         minPts = 9)
umap_dt[, clone := LETTERS[clones$cluster + 1]]


plot = ggplot(umap_dt, aes(x = x, y = y, color = clone)) +
  geom_point(size = 2) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "")

save_and_plot(plot, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/chr11_umap-Clones-", binsize),
              width = 7, height = 7)

# Compare treated vs untreated
# region of interest(roi) is chr11:118,307,205-125,770,541
roi = data.table(chr = "11", start = 118307205, end = 125770541)
roi_bins = foverlaps(roi, bins, which = T)$yid

# Select cells that have any part of the ROI altered (cells of interest; COI)
treated_coi = unlist(pbsapply(colnames(treated[, 4:ncol(treated)]), function(x) {
  count = sum(treated[roi_bins, ..x] != 2)
  if(count > 2) return(x)
}, cl = numthreads))

untreated_coi = unlist(pbsapply(colnames(untreated[, 4:ncol(untreated)]), function(x) {
  count = sum(untreated[roi_bins, ..x] != 2)
  if(count > 2) return(x)
}, cl = numthreads))
  
# Plot
bins[, bin := seq_along(chr)]
bins[, end_cum := cumsum((end - start) + 1)]
bins[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Get bins of -10mb and +10 mb of gRNAs
closest = data.table(chr = "11", start = c(roi$start[1] - 1e7, roi$end[1] + 1e7))
setkey(closest, chr, start)
plot_window = bins[closest, roll = "nearest"]$bin

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

# Subset dt
dt = data.table(cbind(bins, treated[, ..treated_coi], untreated[, ..untreated_coi]))
# dt = cbind(bins, treated, untreated)
# dt = dt[, c(colnames(bins), umap_coi), with = F]
dt = dt[plot_window[1]:plot_window[2], ]

# Order
hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")
dhc = as.dendrogram(hc)

# Rectangular lines
ddata = dendro_data(dhc, type = "rectangle")
annot = data.table(sample = ddata$labels$label)
annot[sample %in% untreated_coi, type := "Untreated"]
annot[sample %in% treated_coi, type := "Treated"]
annot = merge(annot, umap_dt[, c(3, 5)])
annot[, sample := factor(sample, levels = ddata$labels$label)]

# Prepare for heatmap
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 1.8) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), breaks = c(bins[plot_window[1], start_cum],
                                                  bins[roi_bins[1], start_cum],
                                                  bins[roi_bins[length(roi_bins)], start_cum],
                                                  bins[plot_window[2], end_cum]),
                     labels = c("-10 Mb", "gRNA1", "gRNA2", "+10 Mb")) +
  geom_hline(yintercept = bins[roi_bins]$start_cum[c(1, length(roi_bins))], linetype = 2, color = "red") +
  theme(axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

annot_plt1 = ggplot(annot, aes(x = 1, y = sample, fill = type)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d(begin = 0.4) +
  theme_void() +
  theme(legend.position = "left",
        legend.title = element_blank())

annot_plt2 = ggplot(annot, aes(x = 1, y = sample, fill = clone)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_brewer(palette = "Dark2") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank())

combined = plot_grid(annot_plt1, annot_plt2, heatmap, ncol = 3, align = "h", rel_widths = c(0.1, 0.06, 1))
save_and_plot(combined, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/roi_zoomin-", binsize),
              width = 12, height = 5)

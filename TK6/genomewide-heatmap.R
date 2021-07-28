## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting genomeide heatmap of TK6 cells

## Load/install packages
packages = c("data.table", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
numthreads = 32
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

# Combine and sort
total = merge(treated, untreated)
total = total[gtools::mixedorder(total$chr)]
bins = bins[gtools::mixedorder(bins$chr)]

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
dt = cbind(bins, total[, 4:ncol(total)])

# Distance and clustering
hc = hclust(dist(t(dt[, 7:ncol(dt)])), method = "average")

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
  ggrastr::rasterize(geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = .3)) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number", subtitle = paste0("n = ", ncol(dt) -7)) + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())


samples = data.table(sample = colnames(dt[, 7:ncol(dt)]))
samples[, type := ifelse(sample %in% colnames(treated), "Treated", "Untreated")]
samples[, sample := factor(sample, levels = ddata$labels$label)]

annot_plt = ggplot(samples, aes(x = 1, y = sample, fill = type)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d() +
  theme_void() +
  theme(legend.position = "left")

# Plot combined
plt = plot_grid(dendro, annot_plt, heatmap,  align = "h", rel_widths = c(.3, .15, 2), ncol = 3)

save_and_plot(plt, paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/genomewide-", binsize),
              width = 22, height = 14)

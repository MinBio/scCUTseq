## Author: Luuk Harbers
## Date: 2020-11-02
## Script for plotting copy number profiles and genomewide heatmaps

## Load/install packages
packages = c("data.table", "tidyverse", "scales", "pbapply", "ggdendro", "GenomicRanges")
invisible(sapply(packages, require, character.only = T))
invisible(source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R"))

load(snakemake@input[[1]])
data = readRDS(snakemake@input[[2]])

# load("/mnt/AchTeraD/data/ngi/P18158_combined/NZ120/out/dnaobj-psoptim-50000.Rda")
# data = readRDS("/mnt/AchTeraD/data/ngi/P18158_combined/NZ120/out/dnaobj-50000.Rds")

# Select stats
stats = data$stats[[1]]
setDT(stats)

# Select bins
bins = data$binbed[[1]]
setDT(bins)

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

#Plot profile plots
invisible(pblapply(samples, function(sample) {
  dt = cbind(bins, psoptim[, sample], 2^dnalrr[,sample]*psoptim.par[, sample])
  setnames(dt, c("chr", "start", "end", "bin", "end_cum", "start_cum", "cn", "raw"))

  dt[, col := ifelse(cn < 11, as.character(cn), "10+")]
  dt[, col := factor(col, levels = c(as.character(0:10), "10+"))]

  # Sample stats
  stat_string = paste0(sample,
                       " | reads: ", stats[cell == sample,]$reads,
                       " | avg reads/bin: ", as.integer(stats[cell == sample]$mean),
                       " | spikiness: ", round(stats[cell == sample,]$spikiness, 2))

  # save plot
  ggplot(dt, aes(x = bin)) +
    geom_point(aes(y = raw, color = col), size = 0.7) +
    geom_point(aes(y = cn), size = 1) +
    scale_color_manual(values = colors, drop = F) +
    scale_y_continuous(labels=comma_format(accuracy = 1), breaks = pretty_breaks(6)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(y = "Copy Number", x = "", subtitle = stat_string) +
    geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2) +
    geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chr), vjust = -0.5, hjust = "center", inherit.aes = F) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  # Save plot
  ggsave(filename = paste0(snakemake@params[["outdir_profiles"]], snakemake@wildcards[["binsize"]], "/", sample, ".png"),
         width = 14, height = 7, units = "in", dpi = 300)
}, cl = snakemake@threads[[1]]))

# Save empty plots for empty wells
toCreate = snakemake@params[["samples"]][!snakemake@params[["samples"]] %in% samples]
x = data.table()
invisible(lapply(toCreate, function(sample) {
  write.table(x, paste0(snakemake@params[["outdir_profiles"]], snakemake@wildcards[["binsize"]], "/", sample, ".png"),
              col.names = F, row.names = F, quote = F)
}))

# combine data
dt = data.table(cbind(bins, psoptim))

# Make dendrogram
cat("Calculate distances and cluster all samples")
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

# Set sample order
dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

# Plot heatmap 20K+
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
combined = plot_grid(dendro, heatmap,  align = "h", rel_widths = c(0.2, 2), ncol = 2)
save_and_plot(combined, 
              paste0(snakemake@params[["outdir_genomewide"]], "genomewideheatmap_", snakemake@wildcards[["binsize"]]),
              width = 16, height = nrow(ddata$labels) * 0.07, dpi = 900)

# 300K+ samples
# Make dendrogram for 300K+ samples
hqsamples = stats[(reads > 3e5 & mean >= 50 & spikiness < 0.5)]$cell
if(length(hqsamples) > 1) {
  dt = dt[, colnames(dt) %in% c("chr", "start", "end", "bin", "start_cum", "end_cum", hqsamples), with = F]
  
  # Distance and clustering
  cat("Calculate distances and cluster high quality samples")
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
  
  # Set sample order
  dt_melt[, variable := factor(variable, levels = ddata$labels$label)]
  
  # Plot heatmap 20K+
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
  combined = plot_grid(dendro, heatmap,  align = "h", rel_widths = c(.3, 2), ncol = 2)
  save_and_plot(combined, 
                paste0(snakemake@params[["outdir_genomewide"]], "HQ-genomewideheatmap_", snakemake@wildcards[["binsize"]]),
                width = 16, nrow(ddata$labels) * 0.07, dpi = 900)
} else {
    x = data.table()
    write.table(x, paste0(snakemake@params[["outdir_genomewide"]], "HQ-genomewideheatmap_", snakemake@wildcards[["binsize"]], ".png"),
                col.names = F, row.names = F, quote = F)
    }




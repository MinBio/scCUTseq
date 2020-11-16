packages = c("data.table", "ggplot2", "umap", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 16
set.seed = 777

base_path = "/mnt/AchTeraD/data/ngi/P18158_combined/"

rda = list.files(base_path, pattern = "-250000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-250000.Rds", recursive = T, full.names = T)
libs = list.dirs(base_path, recursive = F, full.names = F)


# Set QC thresholds
min_reads = 3e5
max_spikiness = 0.55
min_avgreads = 50
diploid = TRUE

total = pblapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, ..usable_cells]
  
  # Select diploid only
  if(diploid) dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
  
  # Set colnames to include library
  if(ncol(dt) > 0) setnames(dt, paste0(libs[lib], "_", colnames(dt)))
  
  # Return dt
  return(dt)
}, cl = nthreads)

total = do.call(cbind, total)

# Get annotation

# Specify certain alterations to color code on
alterations = data.table(region = c("20q 2-copy gain", "diploid chr13", 
                                    "chr16q deletion", "chr19p deletion", 
                                    "chr19 deletion", "chr11q deletion",
                                    "chr8 amplification"),
                         startbin = c(9369, 7322, 8359, 9066, 9056, 6770, 5366),
                         endbin = c(9445, 7658, 8513, 9138, 9252, 6797, 5454),
                         ploidy = c(4, 2, 1, 1, 1, 1, 3))
alterations[, alt_length := endbin-startbin]

# Get values
clusters = lapply(colnames(total), function(sample) {
  dt = apply(alterations, 1, function(alteration) {
    counts = sum(total[, ..sample][alteration["startbin"]:alteration["endbin"]] == alteration["ploidy"])
    dt = data.table(counts)
    setnames(dt, sample)
    return(dt)
  })
  dt = rbindlist(dt)
  select = dt > alterations$alt_length * 0.75
  annot = alterations[select[, 1]]$region
  

  if("diploid chr13" %in% annot) annot = "diploid chr13"
  else if("chr19 deletion" %in% annot) annot = "chr19 deletion"
  else if("chr16q deletion" %in% annot) annot = "chr16q deletion"
  else if("chr19p deletion" %in% annot) annot = "chr19p deletion"
  else if("chr8 amplification" %in% annot) annot = "chr8 amplification"
  else if("20q 2-copy gain" %in% annot) annot = "20q 2-copy gain"
  #else if("chr19p deletion" %in% annot) annot = "chr19p deletion"
  else if("chr11q deletion" %in% annot) annot = "chr11q deletion"
  else annot = "Major clone"

  #return(data.table(sample = sample, annotation = paste(annot, collapse = ", ")))
  
  return(data.table(sample = sample, annotation = annot))
})

clusters_all = rbindlist(clusters)

# Get annotation
annot = data.table(sample = colnames(total), 
                   library = gsub("_.*", "", colnames(total)),
                   alteration = clusters_all$annotation)


config = umap.defaults
config$n_neighbors = 25
config$min_dist = 0.5

total_umap = umap(t(total), method = "umap-learn", config = config)

umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     group = annot$group,
                     sample = annot$sample,
                     alteration = annot$alteration)

umap_dt[, alteration := factor(alteration, levels = c("Major clone", "20q 2-copy gain", "chr16q deletion",
                                                      "chr19 deletion", "chr8 amplification", "chr11q deletion", 
                                                      "diploid chr13", "chr19p deletion"))]

plot = ggplot(umap_dt, aes(x = x, y = y, color = alteration)) +
  geom_point(size = 2) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Cell type")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-UMAP/TK6-250kb-umap_clustering-alterations",
              height = 6, width = 9)


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

# Sample 5 cells from each alteration
annot = umap_dt[, .SD[sample(.N, min(.N, 5))], by = alteration]
annot[, sample := factor(sample, levels = sample)]
cols = annot$sample

dt = data.table(cbind(bins, total[, ..cols]))

dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
dt_melt[, variable := factor(variable, levels = cols)]

# Plot annotation bar
annot_plt = ggplot(annot, aes(x = 1, y = sample, fill = alteration)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme_void() +
  theme(legend.position = "none")

# Plot heatmap
plt = ggplot(dt_melt) +
  ggrastr::rasterize(geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 5),
                     dpi = 2300) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

combined_plt_noleg = cowplot::plot_grid(annot_plt, plt, ncol = 2, axis = "tb", align = "h", rel_widths = c(0.02, 1))

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = cowplot::plot_grid(combined_plt_noleg, legend, ncol = 1, rel_heights = c(1, 0.08))
# save_and_plot(combined,
#               "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-UMAP/250kb-UMAP-cluster-GenomewideHeatmap",
#               height=7, width=16)
cairo_ps("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-UMAP/250kb-UMAP-cluster-GenomewideHeatmap_rast_cairo_ps.eps",
         onefile = TRUE, width = 16, height = 7, family="Helvetica", pointsize=8, antialias="none", fallback_resolution = 2300)
combined
dev.off()



packages = c("data.table", "ggplot2", "umap", "pbapply", "ggdendro", "ggalt")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 16

base_path = "/mnt/AchTeraD/data/"

rda = list.files(base_path, pattern = "-500000.Rda", recursive = T, full.names = T)
rds = list.files(base_path, pattern = "-500000.Rds", recursive = T, full.names = T)
rda = rda[grepl("BICRO245|BICRO246|BICRO248", rda)]
rds = rds[grepl("BICRO245|BICRO246|BICRO248", rds)]
annot = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/run_barcode_blockid.tsv", header = T)

runs = c("BICRO245", "BICRO246", "BICRO248")

# Set QC thresholds
min_reads = 3e5
max_spikiness = 0.55
min_avgreads = 50

total = pblapply(1:length(rda), function(lib) {
  load(rda[lib])
  stats = readRDS(rds[lib])
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  # Filter out cells from other sequencing
  prost_cells = annot[run == runs[lib], barcode]
  dt = dt[, colnames(dt) %in% prost_cells, with = F]
  
  # Filter out low quality cells
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, colnames(dt) %in% usable_cells, with = F]
  
  # Filter out low quality cells that still manage to sneak through previous filtering
  colmeans = colMeans(dt)
  dt = dt[, colmeans < 5, with = F]
  
  # Set name with block info
  newNames = data.table(run = runs[lib], barcode = colnames(dt))
  newNames = merge(newNames, annot, by = c("run", "barcode"), all.x = T)
  setnames(dt, paste0(newNames$block, "-", newNames$barcode))
  
  # Return dt
  return(dt)
}, cl = nthreads)

# Bind total and get annotations
total = do.call(cbind, total)
blocks = data.table(sample = colnames(total), block = gsub("-.*", "", colnames(total)))

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
colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# Bind dt
dt = cbind(bins, total)

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
  scale_x_continuous(expand = c(0.0, 0.0)) +
  theme_dendro()

# Plot annotation bar
blocks[, sample := factor(sample, levels = ddata$labels$label)]
annot_plt = ggplot(blocks, aes(x = 1, y = sample, fill = block)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme_void() +
  theme(legend.position = "none")

dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]
dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

# Plot heatmap
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 1.2) +
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
                                    align = "h", axis = "tb", rel_widths = c(0.15, 0.01, 1), ncol = 3)

legend = cowplot::get_legend(
  annot_plt + guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

combined = cowplot::plot_grid(combined_noleg, legend, ncol = 1, rel_heights = c(1, 0.05))

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/patient3_genomewideheatmap-BICRO245+246+248",
              height = 20, width = 20)

packages = c("data.table", "ggplot2", "umap", "pbapply", "ggdendro", "ggalt")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

nthreads = 16
set.seed(666)

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

# Make umap
config = umap.defaults
config$n_neighbors = 10  # 20 + 0.5 great
config$min_dist = 0.25

total_umap = umap(t(total), method = "umap-learn", config = config)

umap_dt = data.table(x = total_umap$layout[, 1],
                     y = total_umap$layout[, 2],
                     block = blocks$block,
                     sample = blocks$sample)

plot = ggplot(umap_dt, aes(x = x, y = y, color = block)) +
  geom_point(size = 2) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Block")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/patient3_umap",
              height = 6, width = 9)

# DEFINE CLUSTERS -- MANUALLY
umap_dt[x > 10.1, clusters := "Cluster-1"]
umap_dt[x > 0 & x < 10.1, clusters := "Cluster-2"]
umap_dt[x < 0 & y > -10, clusters := "Cluster-3"]
umap_dt[x > -20 & x < -10 & y > -20 & y < -10, clusters := "Cluster-4"]
umap_dt[x > -20 & x < -15 & y > -21 & y < -17, clusters := "Cluster-5"]
umap_dt[x > -20 & x < -10 & y > -26 & y < -20.5, clusters := "Cluster-6"]
umap_dt[is.na(clusters), clusters := "Cluster-7"]

label_loc = data.table(cluster = as.character(1:7), 
                       x = c(10, 9, -9, -14.5, -18, -20, -15),
                       y = c(-1, -3, -3, -14, -20, -24, -30))

# Plot with circles
plot = ggplot(umap_dt, aes(x = x, y = y, color = block)) +
  geom_point(size = 2) +
  geom_encircle(aes(group = clusters), color = "black", size = 2, expand = 0.01, s_shape = 0.1) +
  geom_text(data = label_loc, aes(x = x, y = y, label = cluster), color = "black", size = 7) +
  scale_color_npg() +
  labs(x  = "UMAP 1", y = "UMAP 2", color = "Block")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/patient3_umap_encircled",
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
colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# sample 8 profiles from each cluster
umap_dt = umap_dt[, .SD[sample(.N, min(.N, 10))], by = clusters]
cols = umap_dt$sample

# Set factors
umap_dt[, sample := factor(sample, levels = sample)]

# Plot annotation bar
annot_plt = ggplot(umap_dt, aes(x = 1, y = sample, fill = clusters)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme_void() +
  theme(legend.position = "left")

dt = cbind(bins, total[, ..cols])

dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

# Set sample order
dt_melt[, variable := factor(variable, levels = umap_dt$sample)]

# Plot heatmap
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 3) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())
  
combined = cowplot::plot_grid(annot_plt,
                              heatmap,
                              align = "h", rel_widths = c(0.08, 1), ncol = 2)

save_and_plot(combined, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/patient3_genomewideheatmap-clusters",
              height = 4, width = 20)

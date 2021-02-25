## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting subclones and getting median CN profile

## Load/install packages
packages = c("data.table", "pbapply", "ggdendro", "tidyr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load in files of clones
clones = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/prostate/P19254-prostate_clones.tsv")
setorder(clones, V2)

# Set threads
nthreads = 32

# Select all libraries with prostate cancer cells)
libraries = data.table(library = c(paste0("MS", 77:81), paste0("NZ", 185:189)),
                       basepath = c(paste0("/mnt/AchTeraD/data/ngi/P19254/MS", 77:81, "/"),
                                    paste0("/mnt/AchTeraD/data/ngi/P19254/NZ", 185:189, "/")))

# Loop through dirs and get copy numbers
total = pblapply(1:nrow(libraries), function(i) {
  load(paste0(libraries$basepath[i], "out/dnaobj-psoptim-500000.Rda"))
  stats = readRDS(paste0(libraries$basepath[i], "out/dnaobj-500000.Rds"))
  bins = stats$binbed[[1]]
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  # Set colnames to include library
  setnames(dt, paste0(libraries$library[i], "_", colnames(dt)))
  
  # Return dt
  return(cbind(bins, dt[, colnames(dt) %in% clones$V1, with = F]))
}, cl = nthreads)

merged = Reduce(function(x, y) merge(x, y, by = c("chr", "start", "end"), all = TRUE, sort = F), total)
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
dt = dt[, c("chr", "start", "end", "bin", "end_cum", "start_cum", clones$V1), with = F]

# Get annotation
clones[, V1 := factor(V1, levels = V1)]
clones[, V2 := factor(V2)]
clones[, Clones := V2]
# clones[V2 == "15", Clones := "A"]
# clones[V2 == "18", Clones := "B"]
# clones[V2 == "19", Clones := "C"]
# clones[V2 == "21", Clones := "D"]
# clones[V2 == "59", Clones := "E"]
# clones[V2 == "176", Clones := "F"]
# clones[V2 == "204", Clones := "G"]
# clones[V2 == "446", Clones := "Diploid"]

annot_plt = ggplot(clones, aes(x = 1, y = V1, fill = Clones)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_viridis_d() +
  theme_void() +
  theme(legend.position = "left")

# Prepare for heatmap
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

# Set sample and chromosome orderorder
dt_melt[, variable := factor(variable, levels = clones$V1)]

# Plot heatmap
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = 1.5) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number", subtitle = paste0("n = ", nrow(clones))) + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

plt = plot_grid(annot_plt, heatmap,  align = "h", rel_widths = c(.12, 2), ncol = 2)

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/prostate/subclones-heatmap",
              width = 20, height = 20)


# Get median cn for each (sub-)clone
dt_clones = merge(dt_melt, clones, by.x = "variable", by.y = "V1")
dt_clones[, value := as.integer(as.character(value))]

# Get median values
clone_median = dt_clones[, .(median_cn = median(value)), by = .(chr, start, end, bin, start_cum, end_cum, V2, Clones)]
clone_median[, median_cn := as.factor(round(median_cn))]
clone_median[, Clones := factor(Clones, levels = c(LETTERS[1:7], "Diploid"))]

heatmap_clones = ggplot(clone_median) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = Clones, color = as.factor(median_cn)), size = 10) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

save_and_plot(heatmap_clones, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/median-cn-subclones-heatmap",
              width = 18, height = 3.3)

# Make wide dt
median_cns = pivot_wider(clone_median, names_from = "V2", values_from = "median_cn")

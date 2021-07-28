## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting subclones and getting median CN profile

## Load/install packages
packages = c("data.table", "pbapply", "ggdendro", "tidyr", "GenomicRanges")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
nthreads = 32

# Select all libraries with breast cancer cells)
# libraries = data.table(library = paste0("NZ", c(249:252, 255, 257)),
#                        basepath = c(paste0("/mnt/AchTeraD/data/BICRO284/NZ", c(249:252, 255, 257), "/")))
libraries = data.table(library = paste0("NZ", 169:174),
                       basepath = c(paste0("/mnt/AchTeraD/data/BICRO277/NZ", 169:174, "/")))
clones = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca_1-clones.tsv")

# Loop through dirs and get copy numbers
total = pblapply(1:nrow(libraries), function(i) {
  rds = readRDS(paste0(libraries$basepath[i], "cnv/500000/cnv.rds"))
  bins = rds$bins
  stats = rds$stats
  dt = data.table(rds$copynumber)
  dt = dt[, colnames(dt) %in% stats[stats$classifier_prediction == "good", sample], with = F]
  # Set colnames to include library
  setnames(dt, paste0(libraries$library[i], "_", colnames(dt)))
  
  # Return dt
  return(cbind(bins, dt))
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
setorder(clones, V2)
clones[, V1 := factor(V1, levels = V1)]
clones[, V2 := as.character(V2)]
# Rename clones
invisible(sapply(1:uniqueN(clones$V2), function(i) {
  clones[, V2 := str_replace(clones$V2, unique(clones$V2)[i], LETTERS[i])]
}))
clones[, Clones := factor(V2)]

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
  scale_fill_hue() +
  theme_void() +
  theme(legend.position = "left")

# Prepare for heatmap
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(as.character(value)) > 10, value := "10+"]
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

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/brca1-subclones-heatmap",
              width = 26, height = 18)


# Get median cn for each (sub-)clone
counts = clones[, .N, Clones]
dt_clones = merge(dt_melt, clones, by.x = "variable", by.y = "V1")
dt_clones = merge(dt_clones, counts, by = "Clones")
dt_clones[, value := as.integer(as.character(gsub("\\+", "", value)))]
#dt_clones[, Clones := paste0(Clones, " (n = ", N, ")")]

# Get median values
clone_median = dt_clones[, .(median_cn = median(value)), by = .(chr, start, end, bin, start_cum, end_cum, V2, N, Clones)]

clone_median[, median_cn := as.factor(round(median_cn))]

heatmap_clones = ggplot(clone_median) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = paste0(Clones, " (n = ", N, ")"), color = as.factor(median_cn)), size = 10) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

save_and_plot(heatmap_clones, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/brca1-median-cn-subclones-heatmap",
              width = 18, height = 6)

# Make wide dt
median_cns = pivot_wider(clone_median, id_cols = c("chr", "start", "end"), names_from = "Clones", values_from = "median_cn")
setDT(median_cns)

# Get shifted dt
shift_dt = as.data.table(data.table::shift(median_cns[, c(1, 4:ncol(median_cns)), with = F], n =-1))
setnames(shift_dt, colnames(median_cns)[c(1, 4:ncol(median_cns))])

# Get breakpoints
indices = unlist(lapply(colnames(shift_dt), function(i) {
  which(median_cns[[i]] != shift_dt[[i]])
}))
indices = unique(sort(c(1, indices + 1)))

medicc_dt = data.table(chr = median_cns[indices, chr],
                       start = median_cns[indices, start],
                       end = median_cns[data.table::shift(indices, -1) , end])
medicc_dt = cbind(medicc_dt, median_cns[indices, 4:ncol(median_cns)])

# Replace NA at end of chr X
medicc_dt[is.na(end), end := median_cns[nrow(median_cns), end]]

# Make long dt for medicc
medicc_dt = pivot_longer(medicc_dt, cols = unique(clone_median$Clones))
setDT(medicc_dt)
setnames(medicc_dt, c("chrom", "start", "end", "sample_id", "cn_a"))

# Add minor allele
medicc_dt[, cn_b := 0]

# Write output for medicc
write.table(medicc_dt[, c(4, 1:3, 5, 6)], "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca1-medicc_dt.tsv", 
            quote = F, col.names = T, row.names = F, sep = "\t")

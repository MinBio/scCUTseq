## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for analysis of Breast cancer cells

## Load/install packages
packages = c("data.table", "pbapply", "copynumber", "purrr", "ggdendro")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
nthreads = 32
removebadbins = 1

# Select all libraries with TK6 cells (control and CAS9)
# libraries = data.table(library = paste0("NZ", 169:174),
#                        basepath = c(paste0("/mnt/AchTeraD/data/ngi/P19254/NZ", 169:174, "/")))
libraries = data.table(library = paste0("NZ", 169:174),
                       basepath = c(paste0("/mnt/AchTeraD/data/ngi/P19254/NZ", 169:174, "/")))
bins = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/CNA-calling/files/hg19/150/variable_500000_150_bwa.bed")
gc = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/CNA-calling/files/hg19/150/GC_variable_500000_150_bwa")

# Loop through dirs and get copy numbers
dt = pblapply(1:nrow(libraries), function(i) {
  counts = fread(paste0(libraries$basepath[i], "out/bincounts-500000.tsv"))

  setnames(counts, paste0(libraries$library[i], "_", colnames(counts)))
  
  return(counts)
}, cl = nthreads)
dt = do.call(cbind, dt)

# Select samples
dt = dt[, colSums(dt) > 2e6, with = F]

# Merge with bins
# merged = cbind(bins[, 1:2], merged)

# Filter out badbins
if(removebadbins==1){
  cat("Calculate outlier bins\n")
  
  pr = rowMeans(sweep(dt, 2, colSums(dt), "/"))
  qt = quantile(pr, probs = 0.995)
  bad = which(pr > qt)
  
  dt = dt[!pr > qt, ]
  bins = bins[!pr > qt, ]
  gc = gc[!pr > qt, ]
}

cat("Normalizing and log2 transformation\n")
lowess.gc = function(x, y) {
  low = lowess(x, log(y), f=0.05);
  z = approx(low$x, low$y, x)
  return(exp(log(y) - z$y))
}

dt = sweep(dt + 1, 2, colMeans(dt + 1), "/")
dt_corrected = pbapply(dt, 2, function(x) {
  lowess.gc(gc[[1]], x + 1 / mean(x + 1))
}, cl = nthreads)

lrr = log2(dt_corrected)

# Mask bins


# Winsorize
normal = cbind(bins[, 1:2], dt_corrected)
wins = winsorize(as.data.frame(normal))

# Segment
joint_segment = multipcf(wins, gamma = 40, Y = as.data.frame(normal))

# Get CN
minPloidy = 1.5
maxPloidy = 6
CNgrid = seq(minPloidy, maxPloidy, by = 0.05)

# Determine Copy Number per sample
cn_res = pblapply(colnames(joint_segment)[6:ncol(joint_segment)], function(sample) {
  # Calculate SoS for CN
  outerRaw = joint_segment[[sample]] %o% CNgrid
  outerRound = round(outerRaw)
  outerDiff = (outerRaw - outerRound) ^ 2
  outerColsums = colSums(outerDiff, na.rm = FALSE, dims = 1)
  CNmult = CNgrid[order(outerColsums)]
  CNerror = round(sort(outerColsums), digits=2)
  CN = CNmult[1]
  
  # Make stats
  list(CN_sos = data.table(CNmult = CNmult, CNerror = CNerror),
       CN = CN,
       profile = round(joint_segment[[sample]] * CN))
}, cl = nthreads)

# Set samplenames
names(cn_res) = colnames(joint_segment)[6:ncol(joint_segment)]

# Get integer copynumbers
cn = lapply(cn_res, function(x) {
  x$profile
})
cn = do.call(cbind, cn)
cn = data.table(cn)
setnames(cn, names(cn_res))

# Plot
segments = joint_segment[, c(1, 3:4)]
setDT(segments)
setnames(segments, c("chr", "start", "end"))

segments[, bin := seq_along(chr)]
segments[, end_cum := cumsum((end - start) + 1)]
segments[, start_cum := c(1, end_cum[1:length(end_cum)-1] + 1)]

# Make chr_bounds
chr_bounds = segments[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end-start)), by = chr]
chr_bounds = chr_bounds %>% 
  mutate(mid = round(min + (max-min) / 2,0),
         end_bp=cumsum(as.numeric(chrlen_bp)), 
         start_bp = end_bp - chrlen_bp, 
         mid_bp = round((chrlen_bp / 2) + start_bp, 0))

#Colors
colors = c("#153570", "#577aba", "#c1c1c1", "#e3b55f", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = c(as.character(0:10), "10+")

# # Select hq samples
# hqsamples = stats[(reads > 3e5 & mean >= 50 & spikiness < 0.5)]$cell

# Make dt
data = data.table(cbind(segments, cn))

# Distance and clustering
hc = hclust(dist(t(data[, 7:ncol(data)])), method = "average")

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
dt_melt = melt(data, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[, value := factor(value)]
dt_melt[as.numeric(value) > 10, value := "10+"]
dt_melt[, value := factor(value, levels = c(as.character(0:10), "10+"))]

# Set sample order
dt_melt[, variable := factor(variable, levels = ddata$labels$label)]

# Calculate required linesize
# linesize = (-1/6) * length(hqsamples) + 17
# linesize = max(linesize, 1.5)

# Plot heatmap
heatmap = ggplot(dt_melt) +
  geom_linerange(aes(ymin = start_cum, ymax = end_cum, x = variable, color = value), size = .5) +
  coord_flip() +
  scale_color_manual(values = colors, drop = F) +
  labs(color = "Copy Number") + 
  scale_y_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid_bp) + 
  geom_hline(data = chr_bounds, aes(yintercept = end_bp), linetype = 1, size = .8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())
# Plot
plot_grid(dendro, heatmap,  align = "h", rel_widths = c(.3, 2), ncol = 2)
plotHeatmap(joint_segment, upper.lim = 4, lower.lim = -4)

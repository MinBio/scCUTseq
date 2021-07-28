## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for analysis of TK6 - CAS9 cells

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Set threads
nthreads = 32

# Select all libraries with TK6 cells (control and CAS9)
libraries = data.table(library = paste0("NZ", c(120:135, 175, 179:183)),
                       basepath = c(paste0("/media/luukharbers/HDD1/P18158_combined/NZ", 120:135, "/"),
                                    paste0("/mnt/AchTeraD/data/ngi/P19254/NZ", c(175, 179:183), "/")),
                       type = c(rep("Cas9", 16), rep("Control", 6)))


# Set thresholds
min_reads = 3e5
max_spikiness = 0.5
min_avgreads = 50
diploid = TRUE

# Loop through dirs and get copy numbers
total = pblapply(1:nrow(libraries), function(i) {
  load(paste0(libraries$basepath[i], "out/dnaobj-psoptim-500000.Rda"))
  stats = readRDS(paste0(libraries$basepath[i], "out/dnaobj-500000.Rds"))
  bins = stats$binbed[[1]]
  stats = stats$stats[[1]]
  setDT(stats)
  dt = data.table(psoptim)
  
  usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
  dt = dt[, ..usable_cells]
  
  # Select diploid only
  if(diploid) dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
  
  # Set colnames to include library
  if(ncol(dt) > 0) setnames(dt, paste0(libraries$library[i], "_", colnames(dt)))
  
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

# Set start and end bin for chr11:118,359,229-125,769,283
delstart = 3233
delend = 3236
plotstart = delstart - 20
plotend = delend + 18

dels = apply(dt[, 7:ncol(dt)], 2, function(x) {
  sum(x[delstart:delend] < 2)
})

dt = dt[plotstart:plotend, c("chr", "start", "end", "bin", "end_cum", "start_cum", names(dels[dels > 2])), with = F]

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
  scale_y_continuous(breaks = c(1690026835, 1700493654, 1702590631, 1711983412),
                     labels = c("-10mb", "gRNA 1", "gRNA2", "Chr11-end")) +
  geom_hline(yintercept = c(1700493654, 1702590631), linetype = 2, color = "red") +
  theme(axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

samples = data.table(sample = colnames(dt[, 7:ncol(dt)]))
samples[, library := gsub("_.*", "", sample)]
annot = merge(samples, libraries, by = "library")
annot[, sample := factor(sample, levels = ddata$labels$label)]

annot_plt = ggplot(annot, aes(x = 1, y = sample, fill = type)) +
  geom_bar(stat = "identity",
           width = 1) +
  scale_fill_npg() +
  theme_void() +
  theme(legend.position = "left")

# Plot combined
plt = plot_grid(annot_plt, heatmap,  align = "h", rel_widths = c(.2, 2), ncol = 3)

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/heatmap-ctrl+cas9_chr11-zoom",
              width = 11, height = 5.3)








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

# Calculate alteration frequency (amp/del) per bin
dt_melt = melt(dt, id.vars = c("chr", "start", "end", "bin", "start_cum", "end_cum"))
dt_melt[value > 2, alteration := "AMP"]
dt_melt[value < 2, alteration := "DEL"]
dt_melt[is.na(alteration), alteration := "Neutral"]

# Annotate with experiment type
dt_melt[, library := gsub("_.*", "", variable)]
dt_melt = merge(dt_melt, libraries, by = "library")

# Get counts
samples = data.table(sample = colnames(dt[, 7:ncol(dt)]))
samples[, library := gsub("_.*", "", sample)]
annot = merge(samples, libraries, by = "library")
counts = table(annot$type)

freqs = dt_melt[, .(amp_freq_cas9 = sum(alteration == "AMP" & type == "Cas9") / counts[1],
                    del_freq_cas9 = sum(alteration == "DEL" & type == "Cas9") / counts[1],
                    amp_freq_cas9 = sum(alteration == "AMP" & type == "Control") * -1 / counts[2],
                    del_freq_cas9 = sum(alteration == "DEL"  & type == "Control") * -1 / counts[2]), 
                by = .(chr, start, end, bin, start_cum, end_cum)]
freqs_m = melt(freqs, id.vars =  c("chr", "start", "end", "bin", "start_cum", "end_cum"))
freqs_m[, cna := gsub("_.*", "", variable)]
freqs_m[, type := gsub(".*_", "", variable)]

# Plot frequencies
plot = ggplot(freqs_m, aes(x = bin, y = value, fill = cna, color = cna, group = type)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid) +
  geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2, size = .8) +
  geom_hline(yintercept = 0, size = .8, linetype = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(y = "Alteration frequency", x = "", fill = "Alteration\ntype")

save_and_plot(plot, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/alteration-frequency_cas9-ctrl",
              width = 20, height = 6)

# Plot frequencies
freqs_m[value > .25, value := .25]
freqs_m[value < -.25 , value := -.25]

plot_zoom = ggplot(freqs_m, aes(x = bin, y = value, fill = cna, color = cna, group = type)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0), labels = chr_bounds$chr, breaks = chr_bounds$mid) +
  geom_vline(data = chr_bounds, aes(xintercept = max), linetype = 2, size = .8) +
  geom_hline(yintercept = 0, size = .8, linetype = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(y = "Alteration frequency", x = "", fill = "Alteration\ntype")

save_and_plot(plot_zoom, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/TK6-deletion/alteration-frequency_cas9-ctrl-yzoom",
              width = 20, height = 6)

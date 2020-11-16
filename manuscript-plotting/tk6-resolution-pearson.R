## Author: Luuk Harbers
## Date: 2020-10-28
## Script for plotting TK6 cells deletion

# Disable scientific notation
options(scipen = 999)

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

blacklist = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/CNA-calling/files/hg19-blacklist.v2_adjusted.bed",
                  col.names = c("chr", "start", "end", "type"))
base_path = "/mnt/AchTeraD/data/ngi/P18158_combined/"
bed_base = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/CNA-calling/files/"
binsizes = c(500000, 250000, 175000, 100000, 50000)
nthreads = 40

total = lapply(binsizes, function(bin) {
  rda = list.files(base_path, pattern = paste0("-", bin, ".Rda"), recursive = T, full.names = T)
  rds = list.files(base_path, pattern = paste0("-", bin, ".Rds"), recursive = T, full.names = T)
  libraries = basename(list.dirs(base_path, recursive = F))
  libraries = libraries[!grepl("fastq", libraries)]
  
  
  # Set thresholds
  min_reads = 3e5
  max_spikiness = 0.55
  min_avgreads = 50
  diploid = TRUE
  
  # Loop through libraries
  data = pblapply(1:length(rda), function(lib) {
    load(rda[lib])
    info = readRDS(rds[lib])

    cn = data.table(psoptim)
    bins = data.table(info$binbed[[1]])
    stats = data.table(info$stats[[1]])

    # Get usable cells and filter dt
    usable_cells = stats[reads > min_reads & spikiness < max_spikiness & mean > min_avgreads, cell]
    dt = cn[, usable_cells, with = F]
    
    # Select majorly diploid only cells
    if(diploid) dt = dt[, apply(dt, 2, function(x) mean(x) > 1.9 & mean(x) < 2.1), with = F]
    
    
    if(ncol(dt) > 0) {
      setnames(dt, paste0(libraries[lib], "_", colnames(dt)))
    }
    
    # return dt
    return(cbind(bins, dt))
  }, cl = nthreads)
  
  # Merge datasets
  total = Reduce(function(...) merge(..., all = TRUE), data)

  # Get correlation
  pairwise = cor(total[, 4:ncol(total)], method = "pearson", use = "complete.obs")
  pairwise[upper.tri(pairwise, diag = T)] = NA
  pairwise = reshape2::melt(pairwise, na.rm = T)
  
  # Return correlation and original data.table
  return(list(correlation = data.table(pearson = pairwise$value, binsize = bin),
              copynumbers = total))
})

cors = rbindlist(lapply(total, function(x) {
  x$correlation
}))

cors[, binsize := paste0(binsize / 1e3, "kb")]
cors[, binsize := factor(binsize, levels = c("500kb", "250kb", "175kb", "100kb", "50kb"))]

obs = cors[, .N, by = binsize]

# Pairwise correlation
ppc = ggplot(cors, aes(x = binsize, y = pearson)) +
  ggrastr::rasterize(geom_violin(aes(fill = binsize, color = binsize)), dpi = 600) +
  geom_boxplot(width = 0.05) + 
  scale_fill_viridis_d(begin = 0.4) + 
  scale_color_viridis_d(begin = 0.4) +
  geom_text(data = obs, aes(y = -.55, label = paste("n =", N))) +
  scale_y_continuous() +
  labs(y = "Pairwise Pearson's Correlation", x = "") +
  theme(legend.position = "none")

save_and_plot(ppc, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-resolutions/TK6_pairwise-pearson-resolutions",
              height = 7, width = 7)

# Select all binsizes but 500kb
ref_cn = total[[1]]$copynumbers
results = lapply(total[2:5], function(x) {
  cn = x$copynumbers
  
  # Apply through samples
  cor = pblapply(colnames(cn[, 4:ncol(cn)]), function(sample) {
    selectCols = c("chr", "start", "end", sample)
    
    # Select sample
    profile = cn[, ..selectCols]
    ref = ref_cn[, ..selectCols]
    
    # Setkeys
    setkey(profile, chr, start, end)
    setkey(ref, chr, start, end)
    
    # Get overlaps
    combined = foverlaps(ref, profile)

    # Collapse smaller binned to large bins
    combined[, feature := paste0(chr, i.start, i.end)]
    collapsed = combined[, lapply(.SD, mean), by = feature, .SDcols = 4]
    
    toCor = cbind(ref, collapsed)
    return(cor(toCor[[4]], toCor[[6]], use = "complete.obs"))
    }, cl = nthreads)
  
  binsize = paste0(x$correlation$binsize[1] / 1e3, "kb")
  
  return(data.table(correlation = unlist(cor), binsize = binsize))
})

total_res = rbindlist(results)

# Plotting
total_res[, binsize := factor(binsize, levels = c("250kb", "175kb", "100kb", "50kb"))]
total_res[, correlation := as.numeric(correlation)]
obs = total_res[, .N, by = binsize]

# Pairwise correlation
mpc = ggplot(total_res, aes(x = binsize, y = correlation)) +
  geom_violin(aes(fill = binsize, color = binsize)) +
  geom_boxplot(width = 0.05) + 
  geom_text(data = obs, aes(y = .6, label = paste("n =", N))) +
  scale_fill_viridis_d(begin = 0.4) + 
  scale_color_viridis_d(begin = 0.4) +
  labs(y = "Matched pearson correlation to reference (500kb)", x = "") +
  theme(legend.position = "none")

save_and_plot(mpc, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/manuscript/TK6-resolutions/TK6_matched-pearson-resolutions_500kbref",
              height = 7, width = 7)

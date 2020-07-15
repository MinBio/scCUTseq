# Author: Luuk Harbers
# Date: 2020-07-01
# Script for pearson correlation calculation of AneuFinder called scCUTseq libraries

# Load/install packages
packages = c("data.table", "pbapply")
for(package in packages){
  if(!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package)
  }
}
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# List the samples and parent directory
sample_name = "BICRO226+227"
parent_dir = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/Aneufinder/"

# Specify binsize, always in format of '1e\\+06', '5e\\+05', etc.
binsize = "5e\\+05"
min_reads = 5e4
# Num_threads
num_threads = 32

# Get libraries and models
libraries = list.dirs(paste0(parent_dir, sample_name), recursive = F, full.names = F)

models = lapply(libraries, function(library) {
  files = list.files(paste0(parent_dir, sample_name, "/", library, "/MODELS/method-dnacopy/"), full.names = T)
  files[grepl("5e\\+05", files)]
})

names(models) = libraries

copynumbers = lapply(models, function(library) {
  pblapply(library, function(cell) {
    load(cell)
    # Get copy number info if > min_reads
    if(model$qualityInfo$total.read.count > min_reads) {
      segments = as.data.table(model$segments)
      setkey(segments, seqnames, start, end)
      bins = as.data.table(model$bins)
      setkey(bins, seqnames, start, end)
      overlaps = foverlaps(bins, segments)[, 9]
      setnames(overlaps, "mean.counts", model$ID)
      return(overlaps)
      }
    }, cl = num_threads)
})

# Bind cells together
copynumbers_library = lapply(copynumbers, function(library) {
  dplyr::bind_cols(library)
})

# Get correlation matrix for heatmap
correlation_dt = lapply(copynumbers_library, function(library) {
  as.data.table(melt(cor(library)))
})

# Get only pairwise correlation numbers
correlation_mat = lapply(copynumbers_library, function(library) {
  mat = cor(library)
  mat[lower.tri(mat, diag = F)]
})
pairwise = reshape2::melt(correlation_mat)

plt1 = ggplot(pairwise, aes(x = L1, y = value)) +
  geom_violin(aes(fill = L1)) +
  geom_boxplot(width = 0.05, outlier.size = 0.9) +
  labs(y = "Pearson's Correlation", x = "") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "")

# Save plot
dir.create(paste0(parent_dir, "../pearson-cor"), recursive = T, showWarnings = F)
save_and_plot(plt1, paste0(parent_dir, "../pearson-cor/", sample_name, "_pairwise-library-pearson-correlation"),
              height = 7, width = 7)

# bind together with library name
correlation_dt = pblapply(1:length(correlation_dt), function(i) {
  correlation_dt[[i]][, Var1 := paste0(names(correlation_dt[i]), "_", Var1)]
  correlation_dt[[i]][, Var2 := paste0(names(correlation_dt[i]), "_", Var2)]
}, cl = num_threads)
correlations = dplyr::bind_rows(correlation_dt)

ggplot(correlations, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c("Pearson's\nCorrelation", option="inferno", begin = 0.4, direction = 1) +
  theme_minimal() +
  theme(axis.text = element_blank())

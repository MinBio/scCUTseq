## Author: Luuk Harbers
## Date: 2020-09-02
## Script for writing unfiltered and filtered readcounts and chr-stops to file
  
## Load/install packages
packages = c("data.table", "rtracklayer")
sapply(packages, require, character.only = T)

binsize = "500kb"
sample = "BICRO235_NZ120"

reads = list.files(paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/hmmcopy/counts/", binsize), pattern = ".wig",
                   full.names = T)
mapp = import.wig(paste0("/mnt/AchTeraD/Documents/references/mappability/hg19-mappability_", binsize, ".wig"))
mapp = as.data.table(mapp)

# # Get chromosome stop indexes
# chr_stops = as.data.table(cumsum(table(mapp$seqnames)))
# chr_stops[, chr := as.character(c(1:22, "X"))]

# Filter out bins that have a mappability < 0.7
filter = which(mapp$score < 0.7)
mapp_filtered = mapp[-filter, ]

# Get filtered chromosome stop indexes
chr_stops_filtered = as.data.table(cumsum(table(mapp_filtered$seqnames)))
chr_stops_filtered[, chr := as.character(c(1:22, "X"))]

# Write chromosome stops
# write.table(chr_stops[, 2:1], 
#             paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/hmmcopy/chromosome-stops-", binsize, ".tsv"),
#             col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(chr_stops_filtered[, 2:1], 
            paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/hmmcopy/chromosome-stops-filtered-", binsize, ".tsv"),
            col.names = F, row.names = F, quote = F, sep = "\t")

# Write filtered read counts to matrix
data_matrix = t(sapply(reads, function(x) {
  data = fread(x, header = F)
  
  data = data[!grepl("fixed", data[[1]]),]
  data = data[-filter, ]
  return(data[[1]])
}))

  write.table(data_matrix, 
              paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/hmmcopy/input_matrix/", sample, "-", binsize, "-inputmatrix.txt"),
              col.names = F, row.names = F, quote = F, sep = " ")

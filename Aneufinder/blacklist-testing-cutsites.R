#Author: Luuk Harbers
#Script to analyse read count distribution over cutsites per bin

#require libraries
require(data.table)
require(AneuFinder)
require(ggplot2)

#load data
cutsites = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/cutsite-distribution/hg19-cutsites.tsv")
bins = fread("/mnt/AchTeraD/Documents/bins/fullbins/hg19_100kbBins_nochr.bed")
bamfile = "/mnt/AchTeraD/data/Aneufinder_refs/merged_HQ_IMR90.bam"
binsize = 100e3 #in bp

#set names and key
setnames(cutsites, c("chr", "start", "end"))
setnames(bins, c("chr", "start", "end"))

setkey(cutsites, "chr", "start", "end")
setkey(bins, "chr", "start", "end")

#get overlaps
overlaps = foverlaps(cutsites, bins, nomatch = NA, mult = "all")

#set key and add 1L, and transform to table
setkey(overlaps, "chr", "start", "end")
overlaps[, id := paste0(chr, ":", start, "-", end)]
overlaps[, value := ifelse(is.na(start), 0L, 1L)]

#get counts
counts = overlaps[, sum(value), id]
counts = cbind(unique(overlaps[, 1:3]), counts)

#add missing bins and give 0L
bins[, id := paste0(chr, ":", start, "-", end)]
bins[, V1 := 0L]
counts = rbind(counts, bins[!bins$id %in% counts$id])

# get readcounts of reference
binned_reads = binReads(file = bamfile, assembly = "hg19", bamindex = bamfile, 
                        chromosomes = c(1:22, "X", "Y"), min.mapq = 30, binsizes = binsize)
binned_reads = as.data.table(binned_reads[[1]])

#get overlaps
setnames(binned_reads, "seqnames", "chr")
setkey(binned_reads, "chr", "start", "end")
setkey(counts, "chr", "start", "end")

data = foverlaps(binned_reads, counts[!is.na(counts$start)])

plt1 = ggplot(data, aes(x = V1, y = counts)) +
  geom_point()

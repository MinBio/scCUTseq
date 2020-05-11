# AUTHOR: Luuk Harbers
# Script to generate the large lists of phase II clinical trials full dataset.

#require packages
require(data.table)
require(naturalsort)
require(pbapply)

#load in data
dots = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/gatk-cnv-output/NZ27_hs/readcounts-tumor/",
                  pattern = "denoisedCR.tsv", full.names = T)
calls = list.files("/mnt/AchTeraD/Documents/Projects/scCUTseq/gatk-cnv-output/NZ27_hs/model_tumor/",
                   pattern = "called.seg", full.names = T)
annot = fread("/mnt/AchTeraD/data/BICRO217/NZ26_NZ27_annot.csv")

#set parameters
outdir = "/mnt/AchTeraD/Documents/Projects/scCUTseq/rds-files/gatk-cnv/"
binsize = "100kb"
sample = "NZ27"

#set cores
cl = 20

# RUN

if(!dir.exists(outdir)) dir.create(outdir)

#get names
barcodes = as.data.table(gsub(".*\\/|.dedup.*", "", calls))
names = merge(barcodes, annot, by.x = "V1", by.y = "BARCODE")

#read in data
calls = lapply(calls, function(x) {
  if(file.info(x)$size == 0) {
    return()
    }
  fread(x, skip = "CONTIG")
})

#name list
names(calls) = paste0(names$V1, "-", names$TYPES)

#read in data
dots = lapply(dots, function(x) {
  if(file.info(x)$size == 0) {
    return()
  }
  fread(x, skip = "CONTIG")
})

#name list
names(dots) = paste0(names$V1, "-", names$TYPES)

#remove NULL entries
calls = Filter(Negate(is.null), calls)
dots = Filter(Negate(is.null), dots)

#subset data and setkeys
calls = pblapply(calls, function(x) {
  setkey(x, CONTIG, START, END)
  x[!grepl("GL|MT|Y", x$CONTIG)]
}, cl = cl)

dots = pblapply(dots, function(x) {
  setkey(x, CONTIG, START, END)
  x[!grepl("GL|MT|Y", CONTIG)]
}, cl = cl)

if(!all.equal(names(dots), names(calls))) {
  calls[names(dots)]
}

#get overlaps and combine into one list
combined = pblapply(1:length(calls), function(x) {
  foverlaps(calls[[x]], dots[[x]])[, c(1:4, 7:9)]
}, cl = cl)
names(combined) = names(calls)

#rm previous lists
rm(calls)
rm(dots)

#set levels and set NA logratios for rows that a ratio of 1 or less bins
combined = lapply(combined, function(x) {
  x[, FEATURE := paste0(CONTIG, ":", START, "-", END)]
  x[NUM_POINTS_COPY_RATIO <= 1, MEAN_LOG2_COPY_RATIO := NA]
})

#save RDS file with all information
saveRDS(combined, paste0(outdir, sample, "_", binsize, ".rds"))

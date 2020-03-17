######## LIBRARIES
require(QDNAseq)
require(foreach)
require(doParallel)

######## PARAMETERS
#set bamfile locations
setwd("/mnt/AchTeraD/data/BICRO203/bamfiles/NZ18/")

#set output directory
output <- "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/MS17/"
binfolder <- "/mnt/AchTeraD/Documents/bins/"
#set number of cores to use
cores = 5

######## RUN

#select files
files <- list.files(
  path = "./",
  pattern = "\\.bam$",
  full.names = T
)

#prepare for parallel computing
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

#slect binsizes (matching in file, thus, "1000kbp", "100kbp", etc.)
binsizes <- c("1000kb", "500kb", "100kb", "50kb", "30kb", "10kb")

for(k in binsizes){
  
  #read bin
  bins <- readRDS(paste0(binfolder, "bin_annotations_", k, ".rds"))
  sampleBinsize <- paste0(k, "/")
  
  #create output directories
  outputdir <- paste0(output, "/", sampleBinsize)
  if(!dir.exists(outputdir)){
    dir.create(outputdir)
    dir.create(paste0(outputdir, "/tsvfiles"), recursive = T)
    dir.create(paste0(outputdir, "/bedfiles"), recursive = T)
  }
  
  
  #run calling
  foreach(x = files) %dopar% {
    
    #load libraries
    library(QDNAseq)
    library(DNAcopy)
    
    readCounts <- binReadCounts(bins,bamfiles=x)  
    readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
    readCountsFiltered <- estimateCorrection(readCountsFiltered)
    copyNumbers <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt",smoothBy = 1L)
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
    copyNumbersCalled <- callBins(copyNumbersSegmented,method = "cutoff")
    #filename<- paste(x, ".pdf", sep="")
    #pdf(filename)
    #plot(copyNumbersCalled)
    #dev.off()
    filename <- paste(x,".bed", sep="")
    exportBins(copyNumbersSegmented, paste0(outputdir, "bedfiles/", filename), format = "bed")  
    #exportBins(copyNumbersCalled, format="vcf")
    filename <- paste(x,".tsv",sep="")
    exportBins(copyNumbersCalled, paste0(outputdir, "tsvfiles/", filename), type = "segments", format="tsv")
  }
}

stopCluster(cl)

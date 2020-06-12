# Author: Luuk Harbers
# Script to calculate coverage of (sc)CUTseq sequence data

packages = c("Rsamtools")
for(package in packages){
  if(!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package)
  }
}

# bamlocation
bamfile = "/mnt/AchTeraD/data/BICRO218/NZ40/bamfiles/ACTGAATC.dedup_q30.bam"

param = ScanBamParam(what = c("pos", "qwidth"))
bam = scanBam(bamfile, param = param)

# Calculate coverage
cov = coverage(IRanges(bam[[1]][["pos"]], width = bam[[1]][["qwidth"]]))

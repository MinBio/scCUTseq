# Aneufinder normalization script

require(BSgenome.Hsapiens.UCSC.hg19)
require(AneuFinder)

#load bamfile thats used as a reference. This should be multiple diploid single cells combined
bamfile = "/mnt/AchTeraD/data/Aneufinder_refs/BICRO231_merged_HQ_IMR90.bam"
#bamfile = "/mnt/AchTeraD/data/BICRO217/NZ26/bamfiles_hs/GACACTCA.dedup_q30.bam"
binsize = 1e5 #in bp

#generate blacklists
bins = binReads(file = bamfile, assembly = "hg19", bamindex = bamfile, 
                chromosomes = c(1:22, "X"), min.mapq = 30, binsizes = binsize)

#plot for visual inspections of read distribution across bins, might have to adjust the quantile parameters for optimal results
lcutoff = quantile(bins$`binsize_1e+05`$counts, 0.1)
ucutoff = quantile(bins$`binsize_1e+05`$counts, 0.99)

plt = plot(bins$`binsize_1e+05`) +
  geom_hline(aes(yintercept = lcutoff), color = "red") +
  geom_hline(aes(yintercept = ucutoff), color = "red")

blacklist = bins$`binsize_1e+05`[bins$`binsize_1e+05`$counts <= lcutoff | bins$`binsize_1e+05`$counts >= ucutoff]
blacklist = reduce(blacklist)

blacklist.file <- "/mnt/AchTeraD/data/Aneufinder_refs/blacklist_merged_BICRO231_HQ_IMR90_1e5"
exportGRanges(blacklist, filename=blacklist.file, header=FALSE,
              chromosome.format='NCBI')

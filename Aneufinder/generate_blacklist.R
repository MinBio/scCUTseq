# Aneufinder normalization script

require(BSgenome.Hsapiens.UCSC.hg19)
require(AneuFinder)

#load bamfile thats used as a reference. This should be multiple diploid single cells combined
bamfile = "/mnt/AchTeraD/data/Aneufinder_refs/merged_HQ_IMR90.bam"
#bamfile = "/mnt/AchTeraD/data/BICRO217/NZ26/bamfiles_hs/GACACTCA.dedup_q30.bam"
binsize = 200e3 #in bp

#generate blacklists
bins = binReads(file = bamfile, assembly = "hg19", bamindex = bamfile, 
                chromosomes = c(1:22, "X"), min.mapq = 30, binsizes = binsize)

#plot for visual inspections of read distribution across bins, might have to adjust the quantile parameters for optimal results
lcutoff = quantile(bins$`binsize_2e+05`$counts, 0.1)
ucutoff = quantile(bins$`binsize_2e+05`$counts, 0.9)

plt = plot(bins$`binsize_2e+05`) +
  geom_hline(aes(yintercept = lcutoff), color = "red") +
  geom_hline(aes(yintercept = ucutoff), color = "red")

blacklist = bins$`binsize_2e+05`[bins$`binsize_2e+05`$counts <= lcutoff | bins$`binsize_2e+05`$counts >= ucutoff]
blacklist = reduce(blacklist)

blacklist.file <- "/mnt/AchTeraD/data/Aneufinder_refs/blacklist_merged_HQ_IMR90_strict"
exportGRanges(blacklist, filename=blacklist.file, header=FALSE,
              chromosome.format='NCBI')

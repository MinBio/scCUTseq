dir = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Scripts/CNA-calling/files/hg19/150/"

bins = list.files(dir, pattern = "^variable", full.names = T)
bnds = list.files(dir, pattern = "^bounds", full.names = T)
gcs = list.files(dir, pattern = "^GC", full.names = T)

lapply(1:length(bins), function(i) {
  bin = fread(bins[i])
  
  bnd = fread(bnds[i])

  gc = fread(gcs[i])
  
  bin[, START := shift(END), by = CHR]
  bed = data.table(bin)
  bed[is.na(START), START := 1]
  bed[, CHR := gsub("chr", "", CHR)]
  
  bin = bin[CHR != "Y",]
  gc = gc[1:nrow(bin),]
  bed = bin[1:nrow(bin) != "Y",]
  bnd = bnd[V1 != "X",]
  
  write.table(bin[, 1:2], bins[i], quote = F, row.names = F, col.names = T, sep = "\t")
  write.table(bed[, c(1, 3, 2)], paste0(bins[i], ".bed"), quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(bnd, bnds[i], quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(gc, gcs[i], quote = F, row.names = F, col.names = F, sep = "\t")
  
})
######## LIBRARIES
require(data.table)
require(naturalsort)
require(ggplot2)
require(parallel)
require(doParallel)

#set sample name
samplename <- "BICRO203-NZ18"

#parameters
binsize = "100kb"

#setwd into directory of bed+tsv files
setwd("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/MS17/100kb/")
#set number of cores to use
cores = 5

#set+create output directory
output <- paste0("/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/profile-plots/", samplename, "/", binsize, "/")
if(!dir.exists(output)){
  dir.create(output, recursive = T)
  dir.create(paste0(output, "/png/"), recursive = T)
  dir.create(paste0(output, "/postscript_eps/"), recursive = T)
  dir.create(paste0(output, "/cairo_eps/"), recursive = T)
}


#load and sort segmented files
files.tsv <- list.files("tsvfiles", pattern = ".tsv$", full.names = TRUE, recursive = FALSE)
files.tsv <- naturalsort(files.tsv)
files.bed <- list.files("bedfiles", pattern = ".bed$", full.names = TRUE, recursive = FALSE)
files.bed <- naturalsort(files.bed)

# #load in possible extra annotation files/selection etc
# annotation <- fread("/home/luukharbers/Documents/Projects/Turin/Data/Turin-samples-annotation.csv", header = F)
# setnames(annotation, c("name", "well", "samplename", "barcode"))
# 
# #grep if necessary
# annotation <- annotation[grepl("MS16|NZ12", annotation$name)]
# files.tsv <- files.tsv[grepl("MS16|NZ12", files.tsv)]
# files.bed <- files.bed[grepl("MS16|NZ12", files.bed)]
# 
# barcodes <- gsub(".dedup.*|.*outdata_|.*III_", "", files.tsv)
# names <- gsub(".*\\/|BICRO.*", "", files.tsv)
# files <- data.table(barcode = barcodes, name = names)
# 
# #merge final annotation file
# annotation <- merge(files, annotation, by = c("name", "barcode"))
# annotation$name <- gsub("16|12", "", annotation$name)

# annotation <- annotation[order(annotation$Barcode)]
# annotation <- fread("E:/SciLife Lab/CUTseq/Revision2/Data/Figure4_selectionNewlabels_corrected.txt")

#############
# #extra selection
# x <- numeric(0)
# 
# for(i in 1:nrow(annotation)){
#   #x[i] <- which(grepl(annotation$`Library ID`[i], files.tsv) & grepl(annotation$Barcode[i], files.tsv)) #non robot
# }
# x <- sort(x)
# 
# files.tsv <- files.tsv[x]
# files.bed <- files.bed[x]
# annotation$runOrder <- as.numeric(gsub("XZ", "", annotation$`Library ID`))
# setorder(annotation, runOrder, Barcode)

#############
#load in full bins 
bins <- fread(paste0("/mnt/AchTeraD/Documents/bins/fullbins/hg19_", binsize, "Bins_nochr.bed"))
setnames(bins, 1:3, c("chromosome", "start", "end"))

#add 1 for consistency between datasets
bins$start <- bins$start + 1

#remove X+Y
bins <- subset(bins, grepl("[0-9]", bins$chromosome))
bins$chromosome <- as.numeric(bins$chromosome)

#setorder
setorder(bins, chromosome, start)

#get x-axis tick locations, switch locations
ticks <- numeric(0)
for(k in 1:22){
  ticks[k] <- nrow(bins[bins$chromosome == k]) / 2
}

#get switch point locations
chrSwitch <- c(0, which(bins$chromosome != dplyr::lag(bins$chromosome)))
labels <- as.character(1:22)

setnames(bins, 1:3, c("chromosome", "start", "end"))
bins$feature <- paste0(bins$chromosome, ":", bins$start, "-", bins$end)

#check if file lengths are equal
if(!length(files.tsv) == length(files.bed)) stop("list of tsv and bed files are not equal in length")

#prepare for parallel computing
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

#loop through files
foreach(i = 1:length(files.tsv)) %dopar%{
  
  require(data.table)
  require(naturalsort)
  require(ggplot2)
  
  #load in data
  data.tsv <- fread(files.tsv[i])
  data.bed <- fread(files.bed[i], header = F, select = 1:5)
  
  setnames(data.tsv, 5, "value")
  setnames(data.bed, c("chromosome", "start", "end", "feature", "value"))
  
  to.add <- subset(bins, !bins$feature %in% data.tsv$feature)
  to.add[, value := NA]
  
  data.tsv <- rbind(data.tsv, to.add)
  data.bed <- rbind(data.bed, to.add)
  
  #make chr numeric
  data.tsv$chromosome <- as.numeric(data.tsv$chromosome)
  data.bed$chromosome <- as.numeric(data.bed$chromosome)
  
  #order
  setorder(data.tsv, chromosome, start)
  setorder(data.bed, chromosome, start)
  
  #change name and factor levels for forced order
  data.tsv$feature <- factor(data.tsv$feature, levels = data.tsv$feature)
  
  data.bed <- data.bed[, c(4, 1:3, 5)]
  data.bed$feature <- factor(data.bed$feature, levels = data.bed$feature)
  
  #data lower or higher than (-)4, set to (-)4
  data.tsv$value[data.tsv$value >= 2] <- 2
  data.tsv$value[data.tsv$value <= -2] <- -2
  
  #set called column for fill color
  data.tsv[data.tsv$value >= log2(2.5/2), color := "red"]
  data.tsv[data.tsv$value <= log2(1.5/2), color := "blue"]
  data.tsv[(data.tsv$value > log2(1.5/2)) & (data.tsv$value < log2(2.5/2)) , color := "black"]
  
  #get sample name
  name <- files.tsv[i]
  name <- gsub(".*/", "", name)
  name <- gsub("\\.bam|\\.tsv", "", name)
  name <- gsub("\\.wgs.*", "", name)
  #for robots
  #name <- paste0(annotation[i, 1], "_", annotation[i, 4])
  
  #for johanBRCA
  #name <- paste0("KI_", annotation[i, 6])
  
  #plot and save segmented line plot
  plt <-
    ggplot(data.tsv, aes(x = feature, y = value, group = 1, color = color)) +
    geom_point(data = data.bed, aes(x = feature, y = value, group = 1), color = "darkgrey", size = 0.1) +
    #geom_point(data = data.bed, aes(x = feature, y = value, group = 1), color = "darkgrey", size = 0.7, alpha = 0.5) +
    geom_point(size = 1.2) +
    geom_vline(xintercept = chrSwitch[2:length(chrSwitch)], linetype = 2, size = 1) +
    ggtitle(name) +
    theme_minimal() +
    #scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
    scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 2)) +
    scale_x_discrete(drop = F, breaks = bins$feature[chrSwitch + ticks], labels = labels) +
    scale_color_identity() +
    theme(text = element_text(family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 10),
          #axis.text = element_blank(),
          line = element_blank(),
          axis.ticks.y = element_line(),
          #axis.ticks.y = element_blank(),
          #axis.line.y = element_blank(),
          axis.line.y = element_line(),
          plot.background = element_rect(fill = "transparent", color = NA))
  
  png(file = paste0(output, "png/", name, ".png"), height = 6, width = 12, units = "in", res = 600)
  print(plt)
  dev.off()
  
  postscript(file = paste0(output, "postscript_eps/", name, "_postscript.eps"), height = 6, width = 12, onefile = TRUE,
             paper="special", pointsize=8, horizontal=TRUE)
  print(plt)
  dev.off()
  
  cairo_ps(file = paste0(output, "cairo_eps/", name, "_cairo.eps"), height = 6, width = 12, onefile = TRUE, pointsize=8)
  print(plt)
  dev.off()
  
  print(paste0("Sample ", i, "/", length(files.tsv), " done."))
}

stopCluster(cl)



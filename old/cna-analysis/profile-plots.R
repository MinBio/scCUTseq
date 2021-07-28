require(data.table)
require(ggplot2)
require(naturalsort)
require(pbapply)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#load in data 
combined = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/rds-files/gatk-cnv/NZ27_100kb.rds")

#outdir
outdir = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/profile-plots/gatk-cnv/"
sample = "NZ27"
binsize = "100kb"

#set cores
cl = 20

# RUN

#get bins for chrom separation in plots
#load in full bins 
bins <- fread(paste0("/mnt/AchTeraD/Documents/bins/fullbins/hg19_", binsize, "Bins_nochr.bed"))
setnames(bins, c("CONTIG", "START", "END"))

#add 1 for consistency between datasets
bins[, START := START + 1]

#remove X+Y
bins <- bins[!grepl("Y", CONTIG)]

#add feature column and set keys and levels
bins[, FEATURE := paste0(CONTIG, ":", START, "-", END)]
bins[, FEATURE := factor(FEATURE, levels = FEATURE)]

#add missing rows in dataset
add = lapply(combined, function(x) {
  bins[!bins$FEATURE %in% x$FEATURE, ]
})

#rbind missing rows
combined = pblapply(1:length(combined), function(x) {
  rbind(combined[[x]], add[[x]], fill = T)
}, cl = cl)
names(combined) = names(add)

#remove 'add' list
rm(add)

#reorder
combined = pblapply(combined, function(x) {
  setkey(x, CONTIG, START, END)
  x[naturalorder(CONTIG)]
}, cl = cl)

combined = pblapply(combined, function(x) {
  x[, FEATURE := factor(FEATURE, levels = FEATURE)]
}, cl = cl)

combined = pblapply(combined, function(x) {
  x[CALL == "+", color := "#E41A1C"]
  x[CALL == "-", color := "#377EB8"]
  x[CALL == "0", color := "black"]
})

#get chromosome switch point locations
chrSwitch <- c(0, which(combined[[1]]$CONTIG != dplyr::lag(combined[[1]]$CONTIG)))

#get x-axis tick locations, switch locations
ticks = sapply(c(as.character(1:22), "X"), function(x) nrow(combined[[1]][CONTIG == x]) / 2)

#generate plots
plots = pblapply(combined, function(x) {
  ggplot() +
    geom_point(data = x, aes(x = FEATURE, y = LOG2_COPY_RATIO), color = "grey", size = 0.5) +
    geom_point(data = x, aes(x = FEATURE, y = MEAN_LOG2_COPY_RATIO, color = color), size = 0.5) +
    geom_vline(xintercept = chrSwitch[2:length(chrSwitch)], linetype = 2, size = 1) +
    scale_y_continuous(limits = c(-2, 4)) + 
    scale_x_discrete(drop = F, breaks = x$FEATURE[chrSwitch + ticks], labels = c(as.character(1:22), "X")) +
    scale_color_identity() +
    theme(text = element_text(family = "Helvetica"),
          legend.position = "",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.y = element_text(size = 10),
          #axis.text = element_blank(),
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          axis.line.y = element_blank())
}, cl = cl)

#save plots
if(!dir.exists(paste0(outdir, sample))){
  dir.create(paste0(outdir, sample))
}

pblapply(names(plots), function(x) {
  ggsave(filename = paste0(outdir, sample, "/", x, ".png"), plot = plots[[x]], dpi = 300, units = "in", height = 10, width = 16)
}, cl = cl)


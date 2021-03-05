## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "gridExtra", "grid", "ggplot2", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#top folder
run_name = "BICRO255"
top_folder = "/mnt/AchTeraD/data/BICRO255/"
save_folder = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/sequencing_QC/BICRO255/"

#list libraries
libraries = list.dirs(top_folder, recursive = F, full.names = F)
libraries = libraries[grepl("NZ|MS", libraries)] #can just change the grepl expression in most cases to match library directories

#get input and output reads after barcode/cutsite extraction
total_reads = data.table(library = character(), variable = numeric(), value = numeric())
sample_reads = data.table(V1 = character(), V2 = numeric(), sample = character(), library = character())

# Loop through libraries
for(library in libraries) {
  
  # Get total reads and demultiplexed reads
  count_file = list.files(paste0(top_folder, library), pattern = "log.txt", full.names = T)
  counts = fread(count_file, sep = ":", header = F)
  
  # Get number of mapped reads and reads mapped within range of cutsite
  cutsitefiles = list.files(paste0(top_folder, library, "/logs"), pattern = "filter", full.names = T)
  cutsitefilter = lapply(cutsitefiles, function(file){
    data = fread(file, sep = ":")
    data[, sample := gsub("-.*|.*\\/", "", file)]
  })
  cutsite_dt = rbindlist(cutsitefilter)
  cutsite_dt[, V1 := ifelse(grepl("pre", V1), "mapped", "mapped in cutsite range")]
  
  # Get deduplicated reads
  dedup_file = list.files(paste0(top_folder, library), pattern = "dedup.tsv", full.names = T)
  dedup = fread(dedup_file)
  mapped_dedup = sum(dedup[[11]])

  

  # Make dt for total reads
  total = data.table(library = library, "total reads" = counts$V2[1], "with barcode" = counts$V2[2],
                     "in cutsite range" = sum(cutsite_dt[V1 == "mapped in cutsite range"]$V2),
                     deduplicated = mapped_dedup )

  total = melt(total, id.vars = "library")
  
  # reads to millions and set factors
  total[, value := value / 1e6]
  
  # Make dt for per sample reads
  dedup_dt = dedup[, c(11, 1)]
  dedup_dt[, V1 := "deduplicated"]
  setnames(dedup_dt, c("V2", "sample", "V1"))
  
  sample = rbind(cutsite_dt, dedup_dt[, c(3, 1, 2)])
  sample[, library := library]
  
  sample[, V1 := factor(V1, levels = c("mapped", "mapped in cutsite range", "deduplicated"))]
  sample[, V2 := V2 / 1e6]
  
  # rbind with total DTs
  total_reads = rbind(total_reads, total)
  sample_reads = rbind(sample_reads, sample)
}

# Set factor levels for total_reads
total_reads[, variable := factor(variable, levels = c("total reads", "with barcode", "in cutsite range", "deduplicated"))]

# # # IF PAIRED RUN THIS
# total_reads[variable == "deduplicated", value := value / 2]
# sample_reads[V1 == "deduplicated", V2 := V2 / 2]

# Plot
plt1 = ggplot(total_reads, aes(x = library, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Reads (millions)",
       x = "") +
  theme(legend.title = element_blank(),
        legend.position = "top")

plt2 = ggplot(sample_reads, aes(x = library, y = V2, color = V1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.7, position = position_jitterdodge()) +
  scale_color_manual(values = brewer.pal(4, "Set1")[2:4], labels = c("mapped", "in cutsite range", "deduplicated")) +
  labs(y = "Reads (millions)",
       x = "")  +
  theme(legend.title = element_blank(),
        legend.position = "top")

# Combine plots
plt = arrangeGrob(plt1, plt2)

# Save plots
if(!dir.exists(save_folder)) dir.create(save_folder, recursive = T)
save_and_plot(grid.draw(plt), paste0(save_folder, "sequence_reads_QC"), width = 12, height = 8)

# Get other statistics for excel sheet
samples_perlib = sapply(libraries, function(x) nrow(sample_reads[library == x & V1 == "mapped"]))
sample_prededup = sapply(libraries, function(x) sum(sample_reads[library == x & V1 == "mapped in cutsite range"]$V2))
sample_postdedup = sapply(libraries, function(x) sum(sample_reads[library == x & V1 == "deduplicated"]$V2))
over300k = sapply(libraries, function(x) nrow(sample_reads[library == x & V1 == "deduplicated" & V2 >= 0.5]))
under300k = sapply(libraries, function(x) nrow(sample_reads[library == x & V1 == "deduplicated" & V2 < 0.3]))
under50k = sapply(libraries, function(x) nrow(sample_reads[library == x & V1 == "deduplicated" & V2 < 0.05]))

# Copy paste this into excel file
sheet_text = cbind(run = run_name, yield = paste0(format(sum(total_reads[variable == "total reads"]$value), digits = 0), " million"), 
                   libraries = libraries, prededup = sample_prededup, 
                   postdedup = sample_postdedup, dup_perc = format((1 - sample_postdedup / sample_prededup ) * 100, digits = 4),
                   over300k = paste0(over300k, "/", samples_perlib), under300k = paste0(under300k, "/", samples_perlib),
                   under50k = paste0(under50k, "/", samples_perlib), success = format(over300k / samples_perlib * 100, digits = 4))
write.table(sheet_text, paste0(top_folder, "library_overview.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")

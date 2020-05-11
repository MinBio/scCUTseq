# Author: Luuk Harbers

# Analysis of sequence reads of multiplexed CUTseq libraries
packages = c("data.table", "gridExtra", "grid", "RColorBrewer")
for(package in packages){
  if(!require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package)
  }
}
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Read in excel files and log file

#top folder
run_name = "BICRO198"
top_folder = "/mnt/AchTeraD/data/BICRO198/"
save_folder = "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/sequencing_QC/BICRO198/"

#list libraries
libraries = list.dirs(top_folder, recursive = F, full.names = F)
libraries = libraries[grepl("NZ", libraries)] #can just change the grepl expression in most cases to match library directories

#get input and output reads after barcode/cutsite extraction
total_reads = data.table(library = character(), input_reads = numeric(), output_reads = numeric())
sample_reads = NULL

# Loop through libraries
for(library in libraries) {
  
  #read in log file with input and output reads after extraction
  extracted = list.files(paste0(top_folder, library), recursive = T, pattern = ".*extracted.log", full.names = T)
  all_reads = fread(extracted, select = 2)
  
  #get input reads
  input_reads = all_reads[[1]][grepl("Input", all_reads[[1]])]
  input_reads = as.numeric(gsub(".*: ", "", input_reads))
  
  #get output reads
  output_reads = all_reads[[1]][grepl("output", all_reads[[1]])]
  output_reads = as.numeric(gsub(".*: ", "", output_reads))
  
  #cbind to data.table
  total_reads = rbind(total_reads, data.table(library = library, input_reads = input_reads, output_reads = output_reads))
  
  #extract per-sample read information for mapping % and duplication
  pre_dedup = fread(paste0(top_folder, library, "/all.tsv"), select = c(1, 7, 9))
  post_dedup = fread(paste0(top_folder, library, "/dedup.tsv"), select = 9)
  
  #make DT of reads info
  reads_dt = data.table(library = library, barcode = pre_dedup[[1]], total = pre_dedup[[2]] + pre_dedup[[3]], 
                            mapped = pre_dedup[[3]], deduplicated = post_dedup[[1]])
  
  #melt and rbind with other samples
  sample_reads = rbind(sample_reads, melt(reads_dt, id.vars = c("library", "barcode")))
}

#plot total reads in combined samples
total_reads = melt(total_reads, id.vars = "library")

#extract combined sample reads
subs_reads = sample_reads[, -"barcode"]
combined_sample = subs_reads[, lapply(.SD, sum), by = .(library, variable)]

total_reads = rbind(total_reads, combined_sample[variable != "total",])
total_reads[, value_m := value / 1e6]

plt1 = ggplot(total_reads, aes(x = library, y = value_m, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1", labels = c("total reads", "reads with barcode and cutsite", "mapped", "deduplicated")) +
  labs(y = "Reads (millions)",
       x = "") +
  theme(legend.title = element_blank(),
        legend.position = "top")

#plot per sample read boxplots
sample_reads[, value_m := value / 1e6]
plt2 = ggplot(sample_reads, aes(x = library, y = value_m, color = variable)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.5, position = position_jitterdodge()) +
  scale_color_manual(values = brewer.pal(4, "Set1")[2:4], labels = c("reads with barcode and cutsite", "mapped", "deduplicated")) +
  labs(y = "Reads (millions)",
       x = "")  +
  theme(legend.title = element_blank(),
        legend.position = "top")

#arrange plots and save
plt = arrangeGrob(plt1, plt2)

if(!dir.exists(save_folder)) dir.create(save_folder, recursive = T)
save_and_plot(grid.draw(plt), paste0(save_folder, "sequence_reads_QC"), width = 12, height = 8)

# Text to copy/paste into the excel sheet 
samples_perlib = sapply(libraries, function(x) nrow(sample_reads[library == x & variable == "total"]))
sample_postdedup = sapply(libraries, function(x) sum(sample_reads[library == x & variable == "deduplicated"]$value))
over300k = sapply(libraries, function(x) nrow(sample_reads[library == x & variable == "deduplicated" & value >= 3e5]))
under300k = sapply(libraries, function(x) nrow(sample_reads[library == x & variable == "deduplicated" & value < 3e5]))
under50k = sapply(libraries, function(x) nrow(sample_reads[library == x & variable == "deduplicated" & value < 5e4]))

# Copy paste this into excel file
sheet_text = cbind(run = run_name, yield = paste0(format(sum(total_reads$input_reads)/1e6, digits = 0), " million"), 
                   libraries = libraries, prededup = total_reads$output_reads, 
                   postdedup = sample_postdedup, dup_perc = format((1 - sample_postdedup / total_reads$output_reads) * 100, digits = 4),
                   over300k = paste0(over300k, "/", samples_perlib), under300k = paste0(under300k, "/", samples_perlib),
                   under50k = paste0(under50k, "/", samples_perlib), succes = format(over300k / samples_perlib * 100, digits = 4))

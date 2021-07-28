## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for plotting subclones and getting median CN profile

## Load/install packages
packages = c("data.table", "pbapply", "uwot", "igraph", "bluster", "GenomicRanges", "dbscan", "tibble", "dplyr", "gtools", "tidyr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Functions
# Run clustering
run_clustering = function(umap_df,
                           k_snn_major = 35,
                           k_snn_minor = 17) {
  
  library(scran)
  # building a snn graph for superclones
  message("Building SNN graph.")
  g_major = scran::buildSNNGraph(umap_df[,c(1:2)], k = k_snn_major, transposed = T)
  
  g_clusters = igraph::membership(igraph::components(g_major))
  g_clusters = paste0("s", g_clusters)
  
  # Clustering
  message("Running hdbscan.")
  subclones = dbscan::hdbscan(umap_df[,c(1:2)],
                               minPts = k_snn_minor)
  umap_df$subclones = paste0("c",subclones$cluster)
  
  # for hdb
  # adding the ones classified as outliers to the closest cluster possible according to euclidean distance
  dist_umap = dist(umap_df[,c(1:2)]) %>% as.matrix() %>% as.data.frame() %>%
    rownames_to_column("cell2") %>%
    gather(key = "cell1",
           value = "dist",
           -cell2) %>%
    dplyr::filter(cell1 != cell2)
  
  dist_min = dist_umap %>%
    right_join(umap_df %>% dplyr::select(cell, subclones), by = c("cell2" = "cell")) %>%
    filter(subclones != "c0") %>%
    group_by(cell1) %>%
    slice_min(dist) %>%
    ungroup()
  
  
  for (i in 1:nrow(umap_df)) {
    
    if(umap_df$subclones[i] == "c0") {
      cellname = rownames(umap_df)[i]
      closest_cell = filter(dist_min, cell1 == rownames(umap_df)[i])$cell2
      closest_cell_cluster = filter(umap_df, cell == closest_cell)$subclones
      umap_df$subclones[i] = closest_cell_cluster
      
    }
    
  }
  
  cl_df = tibble::tibble(
    superclones = g_clusters,
    subclones =  umap_df$subclones,
    cells = rownames(umap_df)
  )
  
  cl_df = cl_df %>%
    arrange(superclones, subclones)
  
  # calculating number of cells in every cluster
  freq_df = janitor::tabyl(umap_df$subclones)
  names(freq_df)[1] = "cluster"
  print(freq_df)
  
  classification = cl_df
  
  message("Done.")
  
  return(classification)
  
}

# Colors
get_colors = function() {
  
  # This function is an unfortunate inheritance from the beginning of this project. Everything uses it
  # so it is was kept.
  
  # subclones
  hues1 = paletteer::paletteer_d("ggsci::default_igv",
                                  n = 35) %>% unclass()
  
  
  colors_vec1 = setNames(hues1, paste0("c", 1:35))
  
  #superclones
  
  colors_g_cl <-
    structure(unclass(paletteer::paletteer_d("yarrr::info", n = 9)),
              names = paste0("s", 1:9))
  colors_g_cl[1] = "#E7A79BFF"
  colors_g_cl[2] = "#90A8C0FF"
  colors_g_cl[6] = "#B4DCF5FF"
  colors_g_cl[7] = "#F2C695FF"
  
  colors_list = list(subclones = colors_vec1,
                      superclones = colors_g_cl)
  
  return(colors_list)
  
}
colors_vector = get_colors()

profiles = readRDS("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/brca2_nondiploid-HQ.rds")

# Set threads
nthreads = 32
umap_neighbours = 4  # brca1 = 40 // brca2 = 4
k_major = 25         # brca1 = 25 // brca2 = 25
k_minor = 5          # brca1 = 11 // brca2 = 5

# Load in profiles

# Run UMAP

set.seed(26)
umap_res = umap(t(profiles[, 4:ncol(profiles)]), 
                metric = "manhattan", 
                min_dist = 0, 
                n_neighbors = umap_neighbours,
                spread = 1,
                n_components = 2,
                n_threads = nthreads)

# Transform into dt
umap_dt = data.table(x = umap_res[, 1],
                     y = umap_res[, 2],
                     cell = as.character(1:nrow(umap_res)))


clusters = run_clustering(umap_dt, k_snn_major = k_major, k_snn_minor = k_minor)
umap_dt = merge(umap_dt, clusters, by.x = "cell", by.y = "cells")

colors_vector_umap = c(colors_vector$subclones, colors_vector$superclones)
plt = ggplot(umap_dt, aes(x = x, y = y)) +
  geom_point(size = 8, alpha = 1, aes(color = superclones)) +
  geom_point(size = 1.5, aes(color = subclones)) +
  scale_color_manual(values = colors_vector_umap) +
  labs(x  = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/brca/navin/brca2-umap_super-subclones",
             height = 7, width = 7)

# Add original names
umap_dt = umap_dt[mixedorder(umap_dt$cell)]
umap_dt[, id := colnames(profiles)[4:ncol(profiles)]]

# Write clones
write.table(umap_dt[, 4:6], "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/brca/navin/brca2clones.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

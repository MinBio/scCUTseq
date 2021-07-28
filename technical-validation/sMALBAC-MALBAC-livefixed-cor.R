## Author: Luuk Harbers
## Date: 2020-10-30
## Script for plotting correlation heatmap

## Load/install packages
packages = c("data.table", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load files
files = c("/mnt/AchTeraD/data/BICRO221/NEBNext/cnv/500000/cnv.rds",
          "/mnt/AchTeraD/data/BICRO229/NEBNext/cnv/500000/cnv.rds")


total = lapply(files, function(i) {
  rds = readRDS(i)
  return(rds$copynumber)
})

dt = do.call(cbind, total)

setnames(dt, c("MALBAC 1:200 - fixed (cell 3)", "MALBAC 1:200 - live (cell 1)",
               "MALBAC 1:200 - fixed (cell 1)", "MALBAC 1:200 - live (cell 2)",
               "MALBAC 1:200 - live (cell 3)", "MALBAC 1:200 - fixed (cell 2)",
               "MALBAC 1:1 - live (cell 4)", "MALBAC 1:1 - fixed (cell 4)",
               "MALBAC 1:1 - fixed (cell 2)", "MALBAC 1:1 - live (cell 1)",
               "MALBAC 1:1 - fixed (cell 3)", "MALBAC 1:1 - fixed (cell 1)",
               "MALBAC 1:1 - live (cell 3)", "MALBAC 1:1 - live (cell 2)"))

res = cor(dt)

res_m = reshape2::melt(res, na.rm = T)
setDT(res_m)

res_m = res_m[grepl("live", Var1) & grepl("fixed", Var2), ]
res_m[, Var1 := factor(Var1, levels = c("MALBAC 1:200 - live (cell 1)", "MALBAC 1:200 - live (cell 2)",
                                        "MALBAC 1:200 - live (cell 3)", "MALBAC 1:1 - live (cell 1)",
                                        "MALBAC 1:1 - live (cell 2)", "MALBAC 1:1 - live (cell 3)",
                                        "MALBAC 1:1 - live (cell 4)"))]
res_m[, Var2 := factor(Var2, levels = rev(c("MALBAC 1:200 - fixed (cell 1)", "MALBAC 1:200 - fixed (cell 2)",
                                            "MALBAC 1:200 - fixed (cell 3)", "MALBAC 1:1 - fixed (cell 1)",
                                            "MALBAC 1:1 - fixed (cell 2)", "MALBAC 1:1 - fixed (cell 3)",
                                            "MALBAC 1:1 - fixed (cell 4)")))]

plt = ggplot(res_m, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_y_discrete("") + 
  scale_x_discrete("", position = "top") +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_viridis("Pearson's\ncorrelation", option="B", begin = 0.75, direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
        axis.line = element_blank())

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/technical-validation/correlation/sMALBAC-MALBAC-livefixed",
              height = 7, width = 8)

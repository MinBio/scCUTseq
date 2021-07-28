dt = data.table(filtering = c("input", "matched_barcode", "mapped in cutsite-range", ">= MAPQ30"),
                "O/N" = c(291776584, 282386370, 241002405, 215362334),
                "4hr" = c(210420654, 204236920, 172917504, 155274355))

dt_melt = melt(dt, id.vars = "filtering")
dt_melt[, filtering := factor(filtering, levels = c("input", "matched_barcode", "mapped in cutsite-range", ">= MAPQ30"))]

plt1 = ggplot(dt_melt, aes(x = variable, y = value, fill = filtering)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  labs(y = "Reads", x = "") +
  theme(legend.position = "none")

dt = data.table(filtering = c("input", "matched_barcode", "mapped in cutsite-range", ">= MAPQ30"),
                "O/N" = c(100, 282386370 / 291776584 * 100, 241002405 / 291776584 * 100, 215362334 / 291776584 * 100),
                "4hr" = c(100, 204236920 / 210420654 * 100, 172917504 / 210420654 * 100, 155274355 / 210420654 * 100))
dt_melt = melt(dt, id.vars = "filtering")
dt_melt[, filtering := factor(filtering, levels = c("input", "matched_barcode", "mapped in cutsite-range", ">= MAPQ30"))]

plt2 = ggplot(dt_melt, aes(x = variable, y = value, fill = filtering)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  labs(y = "Percentage of reads from input", x = "", fill = "") 

plt = cowplot::plot_grid(plt1, plt2, rel_widths = c(.5, .8))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/IVT/4hr-vs-ON-IVT",
              width = 10, height = 5)

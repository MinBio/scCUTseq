## Author: Luuk Harbers
## Date: 2020-08-31
## Script for running HMMcopy

## Load/install packages
packages = c("data.table", "HMMcopy", "argparser", "naturalsort", "ggplot2")
invisible(lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Run HMMcopy on set of files with readcount, gc count and mappability")
parser = add_argument(parser, "--counts", help = "Path to readcounts file")
parser = add_argument(parser, "--gc", help = "Path to gc content file")
parser = add_argument(parser, "--mapp", help = "Path to mappability file")
parser = add_argument(parser, "--e", help = "Parameter used for 'e'. If not provided, will generate a reasonable value")
parser = add_argument(parser, "--strength", help = "Parameter used for 'strength'. If not provided, will generate a reasonable value")
parser = add_argument(parser, "--mu", help = "Parameter used for 'mu'. If not provided, will generate a reasonable value")
parser = add_argument(parser, "--outdir", help = "Path output directory")

# Parse arguments
argv = parse_args(parser)

# Testing purposes
argv = list(counts = "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/hmmcopy-output/BICRO226+227/NZ74+78/readcounts/500kb/reads-ACTGATCG.dedup_q30.bam.wig",
            gc = "/mnt/AchTeraD/Documents/references/gc_content/hg19-gc_500kb.wig",
            mapp = "/mnt/AchTeraD/Documents/references/mappability/hg19-mappability_500kb.wig")

# Make into data.table
reads = wigsToRangedData(argv$counts, argv$gc, argv$mapp)
reads = reads[naturalorder(chr),]
reads[, chr := factor(chr, levels = unique(chr))]

# Normalize reads
#reads_normalized = correctReadcount(reads, mappability = 0.7)
reads_normalized = correctReadcount(reads)

# Set bins that are not ideal to NaN
reads_normalized[is.na(copy) | isFALSE(ideal), copy := NaN]
#reads_normalized[is.na(copy), copy := NaN]

# Set new params
#medians = c(-.9, -0.1, 0.5, 0.8, seq(1.2, 4, .4))
medians = c(-5, seq(-1.4, 4.5, .5)) #SKBR3
names(medians) = as.character(0:12)

params = data.table(strength = 1e30,
                    e = 0.999,
                    mu = medians,
                    lambda = 20,
                    nu = 2.1,
                    kappa = c(25, 670, 670, 100, 25, 25, 25, 10, 5, 5, 5, 5),
                    m = medians,
                    eta = c(5e04, 5e04, 5e05, rep(5e04, 9)),
                    gamma = 3,
                    S = 0.01858295)

# params = HMMsegment(reads_normalized, getparam = T)
# medians = params$mu
reads_segmented = HMMsegment(reads_normalized, param = params)

#reads_segmented$segs
setDT(reads_segmented$segs)

# Set keys
setkey(reads_segmented$segs, chr, start, end)
setkey(reads_normalized, chr, start, end)

reads_total = foverlaps(reads_normalized, reads_segmented$segs)
reads_total[is.nan(copy), median := NA]

# Get uniques
reads_total = unique(reads_total, by = c("chr", "i.start", "i.end"))

reads_total[, feature := factor(paste0(chr, ":", i.start, "-", i.end), levels = paste0(chr, ":", i.start, "-", i.end))]
reads_total[, state := factor(state, levels = as.character(0:11))]

#Colors
colors = c("#496bab", "#9fbdd7", "#c1c1c1", "#e9c47e", "#d6804f", "#b3402e",
           "#821010", "#6a0936", "#ab1964", "#b6519f", "#ad80b9", "#c2a9d1")
names(colors) = as.character(0:11)

ggplot(reads_total, aes(x = feature)) + 
  geom_point(aes(y = copy, color = state), size = 2) +
  scale_color_manual(values = colors, drop = F) +
  geom_point(aes(y = median), color = "black", size = 1) +
  geom_hline(yintercept = medians, linetype = 2) +
  scale_y_continuous(breaks = medians, labels = names(medians), limits = c(medians[1]-1, medians[12]+0.5)) +
  labs(y = "Copy Number", x = "", color = "Copy\nNumber") +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())


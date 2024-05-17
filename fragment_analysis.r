# !/usr/local/bin/Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(dplyr)
  library(reshape2)
  library(BiocGenerics)
  library(GenomicRanges)
  library(stats4)
  library(S4Vectors)
  library(IRanges)
  library(GenomeInfoDb)
})

option_list <- list(
  make_option(c("-s", "--sample"), type = "character", default = NULL, help = "Sample identifier"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output directory for results"),
  make_option(c("--range_start"), type = "integer", default = 100, help = "Fragment length range start"),
  make_option(c("--range_end"), type = "integer", default = 400, help = "Fragment length range end"),
  make_option(c("--step"), type = "integer", default = 5, help = "Step for length calculation"),
  make_option(c("--binsize"), type = "double", default = 5, help = "Size of bins in analysis, in Megabase pairs"),
  make_option(c("--ab"), type = "character", default = NULL, help = "Path to AB.rds file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sample <- opt$sample
workdir <- opt$output
range_start <- opt$range_start
range_end <- opt$range_end
ab_rds <- opt$ab

message(Sys.time(), " - Program Start")

# GC Correction Function
gc.correct <- function(coverage, bias) {
  i <- seq(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE), by = 0.001)
  j <- c(i, max(bias, na.rm = TRUE))
  fit.control <- loess.control(surface = "interpolate", statistics = "approximate", trace.hat = "approximate")
  coverage.trend <- loess(coverage ~ bias, control = fit.control)
  coverage.model <- loess(predict(coverage.trend, j) ~ j)
  coverage.pred <- predict(coverage.model, bias)
  coverage.corrected <- coverage - coverage.pred + median(coverage)
}

AB <- readRDS(ab_rds)

counts <- read.csv(file.path(workdir, "counts_by_width.csv"), header = FALSE)
total_counts <- read.csv(file.path(workdir, "counts.csv"), header = FALSE)[, 1]
bingc <- read.csv(file.path(workdir, "bingc.csv"), header = FALSE)[, 1]
count_gc <- data.frame(total = total_counts)
frag_period_list <- seq(range_start, range_end, by = opt$step)

for (i in 2:length(frag_period_list)) {
  feature_name <- paste("frag", frag_period_list[i - 1], frag_period_list[i] - 1, sep = ".")
  count_gc[, feature_name] <- rowSums(counts[, (frag_period_list[i - 1] - range_start + 3):(frag_period_list[i] - range_start + 2)])
  count_gc[, feature_name] <- gc.correct(count_gc[, feature_name], bingc)
  message(Sys.time(), " - GC correction for ", feature_name)
}

AB$frag.gc <- bingc
for (i in 1:ncol(count_gc)) elementMetadata(AB)[, colnames(count_gc)[i]] <- count_gc[, i]

df <- data.frame(
  id = character(),
  seqnames = character(), arm = character(), nfrag = double(), frag_seq = character(), stringsAsFactors = F
)
for (i in 2:length(colnames(count_gc))) {
  tmp <- df.fr2[, c("id", "seqnames", "arm", colnames(count_gc)[i])]
  colnames(tmp)[4] <- "nfrag"
  tmp.df <- tmp %>%
    group_by(id, seqnames, arm) %>%
    summarize(nfrag = sum(nfrag), .groups = "drop")
  tmp.df$frag_seq <- colnames(count_gc)[i]
  df <- rbind(df, data.frame(tmp.df))
}

df$feature <- paste("FSD", df$seqnames, df$arm, df$frag_seq, sep = ".")

df$nfrag_scale <- scale(df$nfrag)
features.window <- reshape2::dcast(df, id ~ feature, value.var = "nfrag_scale")
features.window <- features.window[df$feature]
features.columns <- c("SampleID", colnames(features.window))
features.window$SampleID <- sample
features.window <- features.window[, features.columns]
# rownames(features.window) = sample
# names(dimnames(features.window)) = c("SampleID", "")
write.csv(features.window, file = file.path(workdir, paste0(sample, ".FSD.csv")), row.names = F, quote = FALSE)
message(Sys.time(), " - done")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(dplyr)
})
opt_list <- list(
  make_option("--indir", type="character", default="."),
  make_option("--perf_out", type="character", default="model_performance_merged.csv"),
  make_option("--auroc_out", type="character", default="auroc_curve_merged.csv")
)
opt <- parse_args(OptionParser(option_list = opt_list))

perf_files <- list.files(opt$indir, pattern="^perf_.*\\.csv$", full.names=TRUE, recursive=TRUE)
auro_files <- list.files(opt$indir, pattern="^auroc_.*\\.csv$", full.names=TRUE, recursive=TRUE)

if (length(perf_files) > 0) {
  perf <- do.call(rbind, lapply(perf_files, read.csv))
  write.table(perf, opt$perf_out, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
}
if (length(auro_files) > 0) {
  au <- do.call(rbind, lapply(auro_files, read.csv))
  write.table(au, opt$auroc_out, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
}


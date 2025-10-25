#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE))
    install.packages("optparse", repos="https://cloud.r-project.org")
  library(optparse)
  if (!requireNamespace("SIAMCAT", quietly = TRUE))
    stop("Please install SIAMCAT first (BiocManager::install('SIAMCAT'))")
  library(SIAMCAT)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--features",     type="character", help="Features CSV"),
  make_option("--meta",         type="character", help="Meta CSV"),
  make_option("--taxa_col",     type="character", default="Genus",
              help="FEATURES taxa/feature ID column (name or 1-based index). Default: Genus"),
  make_option("--meta_id_col",  type="character", default="SampleID",
              help="META sample ID column name. Default: SampleID"),
  make_option("--label_column", type="character",
              help="META label column. Default: Disease_Status"),
  make_option("--case_label",   type="character", default="TB_case",
              help="Positive class label. Default: TB_case"),
  make_option("--covars",       type="character", default=NULL,
              help="Comma-separated list of META columns to keep as covariates (besides label). If omitted, keeps all META cols (except ID)."),
  make_option("--out",          type="character", default="sc_base.rds",
              help="Output RDS path. Default: sc_base.rds")
)))
stopifnot(!is.null(opt$features), !is.null(opt$meta))

msg <- function(...) cat(sprintf(...), "\n")

# -------- Read inputs --------
feature <- read.csv(opt$features, check.names = FALSE, stringsAsFactors = FALSE)
meta    <- read.csv(opt$meta,     check.names = FALSE, stringsAsFactors = FALSE)

# Trim headers (avoids hidden spaces causing mis-matches)
colnames(feature) <- trimws(colnames(feature))
colnames(meta)    <- trimws(colnames(meta))

# -------- FEATURES: use --taxa_col, drop it --------
if (grepl("^[0-9]+$", opt$taxa_col)) {
  id_idx <- as.integer(opt$taxa_col)
  if (id_idx < 1 || id_idx > ncol(feature)) stop("taxa_col index out of range.")
  taxa_name <- colnames(feature)[id_idx]
} else {
  id_idx <- match(opt$taxa_col, colnames(feature))
  if (is.na(id_idx)) stop(sprintf("taxa_col '%s' not found in FEATURES.", opt$taxa_col))
  taxa_name <- opt$taxa_col
}
msg("Using FEATURES taxa column: %s (index %d)", taxa_name, id_idx)

taxa_ids <- trimws(as.character(feature[[id_idx]]))
empty <- which(!nzchar(taxa_ids) | is.na(taxa_ids))
if (length(empty)) taxa_ids[empty] <- paste0("feat_", empty)
rownames(feature) <- make.unique(taxa_ids)

feature <- feature[, -id_idx, drop = FALSE]

# Coerce to numeric & replace NAs
feature[] <- lapply(feature, function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(x)))
feature[is.na(feature)] <- 0

# Collapse duplicate taxa (SIAMCAT requires unique feature IDs)
if (anyDuplicated(rownames(feature))) {
  msg("Found duplicated taxa; collapsing by SUM across samples.")
  feature <- rowsum(as.matrix(feature), group = rownames(feature))
  feature <- as.data.frame(feature, check.names = FALSE)
}

# -------- META: set rownames from --meta_id_col, drop it --------
id_idx_meta <- match(opt$meta_id_col, colnames(meta))
if (is.na(id_idx_meta))
  stop(sprintf("META id column '%s' not found. Available: %s",
               opt$meta_id_col, paste(colnames(meta), collapse=", ")))

# Keep a copy of the sample IDs for later checks
sample_ids <- trimws(as.character(meta[[id_idx_meta]]))
rownames(meta) <- sample_ids
# Drop the ID column; keep frame even if empty
meta <- meta[, setdiff(seq_len(ncol(meta)), id_idx_meta), drop = FALSE]

# -------- Optional: restrict META to selected covariates --------
# Always ensure label column is present
if (!opt$label_column %in% colnames(meta)) {
  stop(sprintf("Label column '%s' not in META (after removing ID). Available: %s",
               opt$label_column, paste(colnames(meta), collapse=", ")))
}

if (!is.null(opt$covars)) {
  requested <- unique(trimws(unlist(strsplit(opt$covars, ","))))
  requested <- requested[nzchar(requested)]
  # Remove label if user listed it; we will add it explicitly to keep order clear
  requested <- setdiff(requested, opt$label_column)

  missing <- setdiff(requested, colnames(meta))
  if (length(missing)) {
    msg("Warning: %d requested covariate(s) not found in META and will be ignored: %s",
        length(missing), paste(missing, collapse=", "))
  }
  keep_covars <- intersect(requested, colnames(meta))
  # Keep label + requested covariates only
  meta <- meta[, c(opt$label_column, keep_covars), drop = FALSE]
  msg("Keeping %d META covariate(s): %s",
      length(keep_covars), if (length(keep_covars)) paste(keep_covars, collapse=", ") else "(none)")
} else {
  # Keep all columns (except ID, which we dropped), i.e., label + any other covariates
  # Ensure label column is first (optional, for readability)
  others <- setdiff(colnames(meta), opt$label_column)
  meta <- meta[, c(opt$label_column, others), drop = FALSE]
  msg("Keeping ALL META columns as covariates (except ID). Count: %d", ncol(meta)-1)
}

# -------- Align samples --------
colnames(feature) <- trimws(colnames(feature))
rownames(meta)    <- trimws(rownames(meta))

common <- intersect(colnames(feature), rownames(meta))
if (!length(common)) stop("No overlapping sample IDs between FEATURES (columns) and META (rownames).")
if (length(common) < ncol(feature))
  msg("Dropping %d feature sample(s) not in META.", ncol(feature) - length(common))
if (length(common) < nrow(meta))
  msg("Dropping %d META row(s) not in FEATURES.", nrow(meta) - length(common))

feature <- feature[, common, drop = FALSE]
meta    <- meta[common, , drop = FALSE]

# -------- Per-sample relative abundance (columns sum to 1) --------
M  <- as.matrix(feature)                 # features x samples
cs <- colSums(M); cs[cs == 0] <- 1
Mrel <- sweep(M, 2, cs, "/")
Mrel[is.na(Mrel)] <- 0

# -------- Labels & SIAMCAT --------
label.meta <- create.label(meta = meta,
                           label = opt$label_column,
                           case  = opt$case_label)

# SIAMCAT expects features x samples
sc.obj <- siamcat(feat = Mrel, label = label.meta, meta = meta)

saveRDS(sc.obj, file = opt$out)
msg("Saved SIAMCAT object to %s", opt$out)


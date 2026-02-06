#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(SIAMCAT)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--model_rds",     type="character", help="Path to trained SIAMCAT RDS (from TRAIN_EVAL)"),
  make_option("--features",      type="character", help="Validation features CSV"),
  make_option("--meta",          type="character", help="Validation metadata CSV"),
  make_option("--label_col",     type="character", default="Disease_Status"),
  make_option("--case_label",    type="character", default="TB_case"),
  make_option("--meta_id_col",   type="character", default="SampleID"),
  make_option("--feat_id_col",   type="character", default="Genus"),
  make_option("--outdir",        type="character", default="."),
  make_option("--prefix",        type="character", default=NULL),
  make_option("--threshold",     type="double",    default=0.5)
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$model_rds) || is.null(opt$features) || is.null(opt$meta)) {
  stop("Missing required arguments: --model_rds, --features, --meta")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
if (is.null(opt$prefix)) {
  opt$prefix <- tools::file_path_sans_ext(basename(opt$model_rds))
}

# infer variant from name (label only)
variant <- if (grepl("(?i)mwmote", tail(strsplit(opt$model_rds, "/")[[1]], 1)) ||
              grepl("(?i)mwmote", opt$prefix)) {
  "mwmote"
} else {
  "siamcat"
}

message(sprintf("[validate] Model RDS: %s | variant guess: %s",
                opt$model_rds, variant))

# ---------- helpers ----------
calculate_metrics <- function(TP, TN, FP, FN){
  total <- TP + TN + FP + FN
  acc   <- if (total > 0) (TP + TN) / total else NA_real_
  sens  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec  <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  prec  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv   <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  f1    <- if (is.finite(prec) && is.finite(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_

  data.frame(TP=TP,TN=TN,FP=FP,FN=FN,
             accuracy=acc,sensitivity=sens,specificity=spec,
             precision=prec,npv=npv,f1=f1, stringsAsFactors=FALSE)
}

# returns metrics ONLY (your original behavior)
calculate_eval_matrix <- function(siamcat_obj, threshold = 0.5) {
  preds <- tryCatch(pred_matrix(siamcat_obj), error=function(e) NULL)
  if (is.null(preds)) preds <- siamcat_obj@pred_matrix
  if (is.null(preds) || nrow(preds) == 0) {
    stop("Prediction matrix is empty in holdout object.")
  }
  labels_vec <- label(siamcat_obj)[["label"]]
  pred_mean <- rowMeans(preds, na.rm = TRUE)
  keep <- is.finite(pred_mean) & is.finite(labels_vec)
  pred_mean <- pred_mean[keep]
  labels_vec <- labels_vec[keep]
  if (!length(pred_mean)) stop("No finite predictions/labels in holdout.")

  predicted_class <- ifelse(pred_mean > threshold, 1L, 0L)
  true_labels     <- ifelse(labels_vec == 1L, 1L, 0L)
  TP <- sum(predicted_class == 1L & true_labels == 1L)
  TN <- sum(predicted_class == 0L & true_labels == 0L)
  FP <- sum(predicted_class == 1L & true_labels == 0L)
  FN <- sum(predicted_class == 0L & true_labels == 1L)
  calculate_metrics(TP, TN, FP, FN)
}

# NEW: export ROC curve points for validation
export_roc_curve <- function(siamcat_obj, out_csv, model_id, variant, threshold = 0.5) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    message("[validate] pROC not installed; skipping ROC export.")
    return(invisible(FALSE))
  }

  preds <- tryCatch(pred_matrix(siamcat_obj), error=function(e) NULL)
  if (is.null(preds)) preds <- siamcat_obj@pred_matrix
  if (is.null(preds) || nrow(preds) == 0) {
    message("[validate] Prediction matrix empty; skipping ROC export.")
    return(invisible(FALSE))
  }

  labels_vec <- label(siamcat_obj)[["label"]]
  pred_mean  <- rowMeans(preds, na.rm = TRUE)

  keep <- is.finite(pred_mean) & is.finite(labels_vec)
  pred_mean  <- pred_mean[keep]
  y          <- ifelse(labels_vec[keep] == 1L, 1L, 0L)

  has_pos  <- any(y == 1L)
  has_neg  <- any(y == 0L)
  nonconst <- (length(unique(pred_mean)) > 1L)

  if (!(has_pos && has_neg && nonconst)) {
    message("[validate] Skipping ROC export (degenerate labels or constant predictions).")
    return(invisible(FALSE))
  }

  roc.obj <- tryCatch(pROC::roc(y, pred_mean, quiet = TRUE), error=function(e) NULL)
  if (is.null(roc.obj)) {
    message("[validate] pROC::roc failed; skipping ROC export.")
    return(invisible(FALSE))
  }

  aucroc_val <- data.frame(
    FPR      = 1 - roc.obj$specificities,
    TPR      = roc.obj$sensitivities,
    model_id = model_id,
    variant  = variant,
    threshold= threshold,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(aucroc_val, out_csv)
  message(sprintf("[validate] Wrote ROC curve points: %s", out_csv))
  invisible(TRUE)
}

choose_col <- function(df, requested, commons) {
  if (requested %in% colnames(df)) return(requested)
  cand <- intersect(commons, colnames(df))
  if (length(cand)) return(cand[1])
  nonnum <- colnames(df)[!sapply(df, is.numeric)]
  if (length(nonnum)) return(nonnum[1])
  colnames(df)[1]
}

#  Reading in the data
VF <- fread(opt$features, check.names = FALSE); VF <- as.data.frame(VF)
VM <- fread(opt$meta,     check.names = FALSE); VM <- as.data.frame(VM)

feat_id_use <- choose_col(VF, opt$feat_id_col,
  c("Genus","Feature","Features","ID","Id","id","Taxon","Species","OTU","ASV","X"))
meta_id_use <- choose_col(VM, opt$meta_id_col,
  c("SampleID","sample_id","Sample","ID","Id","id"))

if (!(opt$label_col %in% colnames(VM))) {
  stop(sprintf("Label column '%s' not found in metadata. Available: %s",
               opt$label_col, paste(colnames(VM), collapse=", ")))
}

rownames(VF) <- VF[[feat_id_use]]; VF[[feat_id_use]] <- NULL
rownames(VM) <- VM[[meta_id_use]]

# relative abundance per sample (columns are samples)
cs <- colSums(VF, na.rm = TRUE); cs[cs == 0] <- 1
VF <- sweep(VF, 2, cs, "/"); VF[is.na(VF)] <- 0

# loading the  model and aligning 
sc.tr <- readRDS(opt$model_rds)

# Validation requires stored models (train.model output)
has_models <- FALSE
try({
  mlist <- sc.tr@model_list
  if (!is.null(mlist) && length(mlist) > 0) has_models <- TRUE
}, silent = TRUE)

if (!has_models) {
  stop(paste0(
    "The loaded SIAMCAT object does not contain trained models (model_list is empty).\n",
    "This is expected for the MWMOTE RDS saved from TRAIN_EVAL (e.g. '*_mwmote.rds'), ",
    "which stores CV predictions but not a reusable fitted model for external validation.\n",
    "Please validate using the plain SIAMCAT model RDS (without '_mwmote' in the name)."
  ))
}

np <- norm_params(sc.tr)
if (is.null(np) || is.null(np[["retained.feat"]]))
  stop("Trained model lacks retained.feat in norm params.")
keep <- np[["retained.feat"]]

miss <- setdiff(keep, rownames(VF))
if (length(miss)) {
  VF <- rbind(VF, matrix(0, nrow=length(miss), ncol=ncol(VF),
                         dimnames=list(miss, colnames(VF))))
}
VF <- VF[keep, , drop=FALSE]

lab  <- create.label(meta=VM, label=opt$label_col, case=opt$case_label)
hold <- siamcat(feat=VF, label=lab, meta=VM)

# holdout prediction 
hold <- tryCatch(
  make.predictions(siamcat = sc.tr, siamcat.holdout = hold, normalize.holdout = TRUE),
  error = function(e) {
    stop(paste0(
      "make.predictions failed on the provided SIAMCAT object.\n",
      "Error was: ", e$message, "\n",
      "If you used an '*_mwmote.rds' object, validate with the plain SIAMCAT model RDS instead."
    ))
  }
)

hold <- evaluate.predictions(hold)

# metrics 
MET <- calculate_eval_matrix(hold, threshold = opt$threshold)
MET[["model_id"]] <- opt$prefix
MET[["variant"]]  <- variant

metrics_file <- file.path(opt$outdir, sprintf("validation_metrics_%s.csv", opt$prefix))
fwrite(MET, metrics_file)

# evaluation plot 
pdf_file <- file.path(opt$outdir, sprintf("validation_evaluation_%s.pdf", opt$prefix))
try(suppressWarnings(model.evaluation.plot(hold, fn.plot = pdf_file)), silent=TRUE)

# ROC curve export 
roc_file <- file.path(opt$outdir, sprintf("validation_roc_%s.csv", opt$prefix))
export_roc_curve(
  siamcat_obj = hold,
  out_csv     = roc_file,
  model_id    = opt$prefix,
  variant     = variant,
  threshold   = opt$threshold
)

message(sprintf("[validate] Done. Metrics: %s", metrics_file))


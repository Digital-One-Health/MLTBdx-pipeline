#!/usr/bin/env Rscript

if (!requireNamespace("ranger", quietly = TRUE))
    install.packages("ranger", repos="https://cloud.r-project.org")
  #library(ranger)
if(!requireNamespace("LiblineaR",quietly=TRUE)) 
   install.packages("LiblineaR", repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(optparse)
  library(SIAMCAT)
  library(data.table)
  library(pROC)
  library(ranger)
  library(LiblineaR)
})

opt_list <- list(
  make_option("--split",            type="character"),
  make_option("--model",            type="character"),
  make_option("--norm",             type="character"),
  make_option("--cutoff",           type="character"),
  make_option("--eval_pdf",         type="character"),
  make_option("--interp_pdf",       type="character"),
  make_option("--auroc_csv",        type="character"),
  make_option("--perf_csv",         type="character"),
  make_option("--model_rds_out",    type="character"),
  make_option("--rf_consens_thresh", type="double"),
  make_option("--thresh",           type="double", default=0.5)
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$split), !is.null(opt$model))

msg <- function(...) cat(sprintf("[train_eval] %s\n", sprintf(...)))

# ---------- metric helpers ----------
calculate_metrics <- function(TP, TN, FP, FN){
  total <- TP + TN + FP + FN
  acc   <- if (total > 0) (TP + TN) / total else NA_real_
  sens  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec  <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  prec  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  f1    <- if (is.finite(prec) && is.finite(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_
  mcc_denom <- suppressWarnings(sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
                                     as.numeric(TN + FP) * as.numeric(TN + FN)))
  mcc   <- if (is.na(mcc_denom) || mcc_denom == 0) NA_real_ else ((TP * TN) - (FP * FN)) / mcc_denom
  N  <- total
  Po <- if (N > 0) (TP + TN) / N else NA_real_
  Pe <- if (N > 0) (((TP + FP) * (TP + FN) + (FN + TN) * (FP + TN)) / (N^2)) else NA_real_
  kappa <- if (is.finite(Pe) && (1 - Pe) != 0) (Po - Pe) / (1 - Pe) else NA_real_

  data.frame(
    tp=TP, tn=TN, fp=FP, fn=FN,
    acc=acc, sens=sens, spec=spec,
    precision=prec, f1=f1, kappa=kappa, mcc=mcc,
    stringsAsFactors = FALSE
  )
}

safe_pred_matrix <- function(sc){
  pm <- tryCatch(pred_matrix(sc), error=function(e) NULL)
  if (is.null(pm)) pm <- sc@pred_matrix
  pm
}

# ---------- load split object ----------
msg("Loading split RDS: %s", opt$split)
sc0 <- readRDS(opt$split)

# sanity: retained features must exist
np <- tryCatch(norm_params(sc0), error=function(e) NULL)
rf <- if (!is.null(np)) np[["retained.feat"]] else NULL
if (is.null(rf) || length(rf) == 0) {
  stop(sprintf("No retained features after filtering for norm=%s, cutoff=%s. Lower your cutoff.",
               opt$norm, opt$cutoff))
}

# sanity: labels must have both classes
lab_vec <- label(sc0)[["label"]]
tbl <- table(lab_vec, useNA="no")
msg("Label table: %s", paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse=", "))
if (length(unique(lab_vec)) < 2) {
  stop("Only one label class present in training split; cannot evaluate.")
}

# ---------- train, predict, evaluate ----------
t0 <- proc.time()[3]
sc1 <- train.model(sc0, method = opt$model)
sc1 <- make.predictions(sc1)
sc1 <- evaluate.predictions(sc1)
t1 <- proc.time()[3]

# get probabilities and labels, drop any NA rows
pm <- safe_pred_matrix(sc1)
if (is.null(pm) || nrow(pm) == 0) stop("Prediction matrix is empty after training.")
pred_mean <- rowMeans(pm, na.rm = TRUE)
lbl <- label(sc1)[["label"]]
keep <- is.finite(pred_mean) & is.finite(lbl)
pred_mean <- pred_mean[keep]
lbl       <- lbl[keep]
if (length(pred_mean) == 0) stop("All predictions or labels are NA after filtering.")

# threshold classification
thr <- opt$thresh
yhat <- ifelse(pred_mean > thr, 1L, 0L)
y    <- ifelse(lbl == 1L, 1L, 0L)
TP <- sum(yhat==1L & y==1L); TN <- sum(yhat==0L & y==0L)
FP <- sum(yhat==1L & y==0L); FN <- sum(yhat==0L & y==1L)

met <- calculate_metrics(TP, TN, FP, FN)
met$norm   <- opt$norm
met$model  <- opt$model
met$cutoff <- opt$cutoff
met$model_id <- sprintf("%s_%s_%s", opt$norm, opt$model, opt$cutoff)
met$train_time_sec <- t1 - t0

# ---------- write perf row ----------
perf_out <- as.data.table(met)
setcolorder(perf_out, c("norm","model","cutoff","model_id",
                        "tp","tn","fp","fn","acc","sens","spec","precision",
                        "f1","kappa","mcc","train_time_sec"))
fwrite(perf_out, opt$perf_csv, append = file.exists(opt$perf_csv))

# ---------- AUROC curve (if both classes and non-constant scores) ----------
has_pos <- any(y==1L); has_neg <- any(y==0L)
nonconst <- (length(unique(pred_mean)) > 1L)
if (has_pos && has_neg && nonconst) {
  roc.obj <- tryCatch(roc(y, pred_mean, quiet=TRUE), error=function(e) NULL)
  if (!is.null(roc.obj)) {
    aucroc <- data.frame(
      FPR   = 1 - roc.obj$specificities,
      TPR   = roc.obj$sensitivities,
      model = opt$model,
      filter = opt$cutoff,
      stringsAsFactors = FALSE
    )
    fwrite(aucroc, opt$auroc_csv, append = file.exists(opt$auroc_csv))
  } else {
    msg("pROC::roc failed; skipping AUROC curve write.")
  }
} else {
  msg("Skipping AUROC curve (classes or scores degenerate).")
}

# ---------- plots (best-effort) ----------
try(suppressWarnings(model.evaluation.plot(sc1, fn.plot = opt$eval_pdf)), silent=TRUE)
if (opt$model == 'randomForest'){
#try(suppressWarnings(model.interpretation.plot(sc1, fn.plot = opt$interp_pdf,
#                                               consens.thres = 0.5, limits = c(-3,3),
#                                               heatmap.type='zscore')), silent=FALSE)
tryCatch({
    model.interpretation.plot(
      sc1, fn.plot = opt$interp_pdf,
      consens.thres = opt$rf_consens_thresh,
      limits = c(-3, 3), heatmap.type = "zscore"
    )
  }, error = function(e) {
    message(sprintf("[interp] consens.thres=%s failed (%s). Retrying with 0.", opt$rf_consens_thresh, e$message))
    try(
      suppressWarnings(model.interpretation.plot(
        sc1, fn.plot = opt$interp_pdf,
        consens.thres = 0,
        limits = c(-3, 3), heatmap.type = "zscore"
      )),
      silent = TRUE
    )
  })


}else {
    try(suppressWarnings(model.interpretation.plot(sc1, fn.plot = opt$interp_pdf,
                                               consens.thres = 0.5 , limits = c(-3,3),
                                               heatmap.type='zscore')), silent=FALSE)
}


# ---------- save trained model ----------
saveRDS(sc1, opt$model_rds_out)
msg("Done: %s", met$model_id)


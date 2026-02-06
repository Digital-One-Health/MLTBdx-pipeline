#!/usr/bin/env Rscript

## installing dependencies
if (!requireNamespace("ranger", quietly = TRUE))
    install.packages("ranger", repos="https://cloud.r-project.org")
if (!requireNamespace("LiblineaR", quietly = TRUE))
    install.packages("LiblineaR", repos="https://cloud.r-project.org")
if (!requireNamespace(c("imbalance","xgboost"), quietly= TRUE))
    install.packages(c("imbalance","xgboost"),quietly=TRUE)
suppressPackageStartupMessages({
  library(optparse)
  library(SIAMCAT)
  library(data.table)
  library(pROC)
  library(ranger)
  library(LiblineaR)
  library(imbalance)
  library(glmnet)
  library(xgboost)
})

## command-line options
opt_list <- list(
  make_option("--split",             type="character"),
  make_option("--model",             type="character"),
  make_option("--norm",              type="character"),
  make_option("--cutoff",            type="character"),
  make_option("--eval_pdf",          type="character"),
  make_option("--interp_pdf",        type="character"),
  make_option("--auroc_csv",         type="character"),
  make_option("--perf_csv",          type="character"),
  make_option("--model_rds_out",     type="character"),
  make_option("--rf_consens_thresh", type="double"),
  make_option("--thresh",            type="double", default=0.5),
  make_option("--seed",              type="double",default=42.0)
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$split), !is.null(opt$model))

msg <- function(...) cat(sprintf("[train_eval] %s\n", sprintf(...)))

## metric computation functions
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

## MWMOTE and independent ML functions
siamcat_mwmote_with_custom_ml <- function(sc.obj,
                                          opt,
                                          minority.class = "case",
                                          majority.class = "control",
                                          kMinority = 5,
                                          seed =opt$seed ) {

  if (is.null(sc.obj)) {
    stop("sc.obj is NULL. Please pass a valid SIAMCAT object.")
  }
  if (!methods::is(sc.obj, "siamcat")) {
    stop("sc.obj is not a SIAMCAT object. Did you pass the result of siamcat()?")
  }

  set.seed(seed)

  # Get normalized features and labels
  X_norm <- get.norm_feat.matrix(sc.obj)   
  sample_ids <- colnames(X_norm)
  X <- t(X_norm)                          

  lab <- label(sc.obj)
  y_si <- ifelse(lab$label == 1, minority.class, majority.class)
  y_si <- factor(y_si, levels = c(minority.class, majority.class))

  # CV split
  ds <- data_split(sc.obj)
  num.folds    <- ds$num.folds
  num.resample <- ds$num.resample
  train_all    <- ds$training.folds
  test_all     <- ds$test.folds

  n_samples <- nrow(X)
  pred_mat <- matrix(NA_real_, nrow = n_samples, ncol = num.resample)
  rownames(pred_mat) <- rownames(X)
  colnames(pred_mat) <- paste0("CV_rep", seq_len(num.resample))

  plot_mwmote <- !is.null(opt$plot_mwmote) && isTRUE(opt$plot_mwmote)
  method <- opt$model   

  # Resample Loop
  for (r in seq_len(num.resample)) {

    message("Resample ", r, " / ", num.resample)

    train_folds <- train_all[[r]]
    test_folds  <- test_all[[r]]

    # Fold Loop
    for (k in seq_len(num.folds)) {

      message("  Fold ", k, " / ", num.folds)

      # fold indices 
      train_idx <- train_folds[[k]]
      test_idx  <- test_folds[[k]]

      X_train <- X[train_idx, , drop = FALSE]
      y_train <- y_si[train_idx]

      # Training dataframe for MWMOTE 
      train_df_mwm <- as.data.frame(X_train)
      train_df_mwm$class <- y_train

      n_minority <- sum(train_df_mwm$class == minority.class)
      n_majority <- sum(train_df_mwm$class == majority.class)
      n_syn      <- max(0, n_majority - n_minority)

      # Before MWMOTE plot
      if (plot_mwmote) {
        Xp <- X_train
        sds <- apply(Xp, 2, sd)
        Xp  <- Xp[, sds > 0, drop = FALSE]

        if (ncol(Xp) >= 2) {
          pc <- prcomp(Xp, scale. = TRUE)
          pc_df <- data.frame(pc$x[, 1:2, drop = FALSE], class = y_train)

          par(mfrow = c(1, 2))
          plot(pc_df$PC1, pc_df$PC2,
               col  = ifelse(pc_df$class == minority.class, "red", "blue"),
               pch  = 19,
               main = paste("Before MWMOTE (Rep", r, "Fold", k, ")"),
               xlab = "PC1", ylab = "PC2")
        } else {
          message("    Skipping BEFORE-MWMOTE PCA in rep ", r,
                  ", fold ", k,
                  " (less than 2 non-constant features).")
        }
      }

      # APPLY MWMOTE
      if (n_syn > 0 && n_minority > 0) {
        syn <- mwmote(
          dataset      = train_df_mwm,
          numInstances = n_syn,
          kNoisy       = 5,
          kMajority    = 3,
          kMinority    = kMinority,
          threshold    = 5,
          cmax         = 2,
          cclustering  = 3,
          classAttr    = "class"
        )

        syn <- as.data.frame(syn)

        if (!"class" %in% colnames(syn)) {
          colnames(syn)[ncol(syn)] <- "class"
        }

        if (ncol(syn) != ncol(train_df_mwm)) {
          stop("mwmote returned ", ncol(syn),
               " columns but training data has ",
               ncol(train_df_mwm),
               ". Please inspect colnames(train_df_mwm) and colnames(syn).")
        }

        colnames(syn) <- colnames(train_df_mwm)
        train_bal <- rbind(train_df_mwm, syn)

      } else {
        train_bal <- train_df_mwm
      }

      # After MWMOTE plot
      if (plot_mwmote) {
        Xp2 <- as.matrix(subset(train_bal, select = -class))
        sds2 <- apply(Xp2, 2, sd)
        Xp2  <- Xp2[, sds2 > 0, drop = FALSE]

        if (ncol(Xp2) >= 2) {
          pc2 <- prcomp(Xp2, scale. = TRUE)
          pc2_df <- data.frame(pc2$x[, 1:2, drop = FALSE], class = train_bal$class)

          plot(pc2_df$PC1, pc2_df$PC2,
               col  = ifelse(pc2_df$class == minority.class, "red", "blue"),
               pch  = 19,
               main = paste("After MWMOTE (Rep", r, "Fold", k, ")"),
               xlab = "PC1", ylab = "PC2")
        } else {
          message("    Skipping AFTER-MWMOTE PCA in rep ", r,
                  ", fold ", k,
                  " (less than 2 non-constant features).")
        }
      }

      # Prepare ML inputs 
      X_bal <- as.matrix(subset(train_bal, select = -class))
      y_bal01 <- ifelse(train_bal$class == minority.class, 1, 0)
      y_bal_fac <- factor(y_bal01, levels = c(0, 1))

      # Model Selection 
      if (method == "randomForest") {

        if (!requireNamespace("ranger", quietly = TRUE)) {
          stop("Please install 'ranger'")
        }

        mtry_val   <- if (!is.null(opt$mtry)) opt$mtry else floor(sqrt(ncol(X_bal)))
        ntree_val  <- if (!is.null(opt$num_trees)) opt$num_trees else 500

        df_rf <- data.frame(label = y_bal_fac, X_bal, check.names = FALSE)

        fit <- ranger::ranger(
          dependent.variable.name = "label",
          data = df_rf,
          probability = TRUE,
          mtry = mtry_val,
          num.trees = ntree_val,
          importance = "permutation",
          classification = TRUE
        )

        pred_fun <- function(object, newdata) {
          pr <- predict(
            object,
            data = data.frame(newdata, check.names = FALSE)
          )$predictions
          as.numeric(pr[, "1"])   # P(class==1) = minority
        }

      } else if (method %in% c("ridge", "lasso", "enet")) {

        if (!requireNamespace("glmnet", quietly = TRUE)) {
          stop("Please install 'glmnet'")
        }

        alpha <- if (method == "lasso") {
          1.0
        } else if (method == "ridge") {
          0.0
        } else {
          if (!is.null(opt$alpha)) opt$alpha else 0.5
        }

        fit <- glmnet::cv.glmnet(
          x      = X_bal,
          y      = y_bal01,
          family = "binomial",
          alpha  = alpha
        )

        pred_fun <- function(object, newdata) {
          as.numeric(
            predict(
              object,
              newx = as.matrix(newdata),
              s = "lambda.min",
              type = "response"
            )
          )
        }

      } else if (method %in% c("ridge_ll", "lasso_ll")) {

        if (!requireNamespace("LiblineaR", quietly = TRUE)) {
          stop("Please install 'LiblineaR'")
        }

        type_ll <- if (method == "lasso_ll") 6 else 0

        fit <- LiblineaR::LiblineaR(
          data   = X_bal,
          target = as.integer(y_bal01),
          type   = type_ll,
          cost   = 1,
          bias   = 1,
          cross  = 0,
          verbose = FALSE
        )

        pred_fun <- function(object, newdata) {
          pr <- predict(object, as.matrix(newdata), proba = TRUE)
          as.numeric(pr$probabilities[, 2])  # P(class==1)
        }

      } else if (method == "xgboost") {

        if (!requireNamespace("xgboost", quietly = TRUE)) {
          stop("Please install 'xgboost'")
        }

        X_xgb <- as.matrix(subset(train_bal, select = -class))
        y_xgb <- as.numeric(train_bal$class == minority.class)

        dtrain <- xgboost::xgb.DMatrix(
          data  = X_xgb,
          label = y_xgb
        )

        nrounds <- if (!is.null(opt$nrounds)) opt$nrounds else 200

        fit <- xgboost::xgboost(
          data      = dtrain,
          objective = "binary:logistic",
          nrounds   = nrounds,
          verbose   = 0
        )

        pred_fun <- function(object, newdata) {
          as.numeric(
            predict(object, newdata = as.matrix(newdata))
          )
        }

      } else {
        stop("Unknown method for MWMOTE pipeline: ", method)
      }

      # Predict on test fold for this resample 
      X_test <- X[test_idx, , drop = FALSE]
      p_case <- pred_fun(fit, X_test)

      pred_mat[test_idx, r] <- p_case
    } # end fold loop
  } # end resample loop

  ##  Push predictions back to SIAMCAT
  pred_mat <- pred_mat[sample_ids, , drop = FALSE]
  pred_matrix(sc.obj) <- pred_mat

  return(sc.obj)
}

# load split object 
msg("Loading split RDS: %s", opt$split)
sc0 <- readRDS(opt$split)

# Ensuring retained features must exist
np <- tryCatch(norm_params(sc0), error=function(e) NULL)
rf <- if (!is.null(np)) np[["retained.feat"]] else NULL
if (is.null(rf) || length(rf) == 0) {
  stop(sprintf("No retained features after filtering for norm=%s, cutoff=%s. Lower your cutoff.",
               opt$norm, opt$cutoff))
}

#Cross checking that the labels must have both classes
lab_vec <- label(sc0)[["label"]]
tbl <- table(lab_vec, useNA="no")
msg("Label table: %s", paste(sprintf("%s=%d", names(tbl), as.integer(tbl)), collapse=", "))
if (length(unique(lab_vec)) < 2) {
  stop("Only one label class present in training split; cannot evaluate.")
}

## SIAMCAT train, predict, evaluate
msg("Running plain SIAMCAT with method: %s", opt$model)
t0 <- proc.time()[3]
sc1 <- train.model(sc0, method = opt$model)
sc1 <- make.predictions(sc1)
sc1 <- evaluate.predictions(sc1)
t1 <- proc.time()[3]

pm <- safe_pred_matrix(sc1)
if (is.null(pm) || nrow(pm) == 0) stop("Prediction matrix is empty after training (plain).")
pred_mean <- rowMeans(pm, na.rm = TRUE)
lbl <- label(sc1)[["label"]]
keep <- is.finite(pred_mean) & is.finite(lbl)
pred_mean <- pred_mean[keep]
lbl       <- lbl[keep]
if (length(pred_mean) == 0) stop("All predictions or labels are NA after filtering (plain).")

thr <- opt$thresh
yhat <- ifelse(pred_mean > thr, 1L, 0L)
y    <- ifelse(lbl == 1L, 1L, 0L)
TP <- sum(yhat==1L & y==1L); TN <- sum(yhat==0L & y==0L)
FP <- sum(yhat==1L & y==0L); FN <- sum(yhat==0L & y==1L)

met_plain <- calculate_metrics(TP, TN, FP, FN)
met_plain$norm   <- opt$norm
met_plain$model  <- opt$model
met_plain$cutoff <- opt$cutoff
met_plain$model_id <- sprintf("%s_%s_%s", opt$norm, opt$model, opt$cutoff)
met_plain$train_time_sec <- t1 - t0
met_plain$variant <- "siamcat"

# plots for plain SIAMCAT (with interpretation)
try(suppressWarnings(model.evaluation.plot(sc1, fn.plot = opt$eval_pdf)), silent=TRUE)
if (opt$model == 'randomForest') {
  tryCatch({
    model.interpretation.plot(
      sc1, fn.plot = opt$interp_pdf,
      consens.thres = opt$rf_consens_thresh,
      limits = c(-3, 3), heatmap.type = "zscore"
    )
  }, error = function(e) {
    message(sprintf("[interp] consens.thres=%s failed (%s). Retrying with 0.",
                    opt$rf_consens_thresh, e$message))
    try(
      suppressWarnings(model.interpretation.plot(
        sc1, fn.plot = opt$interp_pdf,
        consens.thres = 0,
        limits = c(-3, 3), heatmap.type = "zscore"
      )),
      silent = TRUE
    )
  })
} else {
  tryCatch({
    model.interpretation.plot(
      sc1, fn.plot = opt$interp_pdf,
      consens.thres = 0.5,
      limits = c(-3, 3), heatmap.type = "zscore"
    )
  }, error = function(e) {
    message(sprintf("[interp] consens.thres=0.5 failed (%s). Retrying with 0.", e$message))
    try(
      suppressWarnings(model.interpretation.plot(
        sc1, fn.plot = opt$interp_pdf,
        consens.thres = 0,
        limits = c(-3, 3), heatmap.type = "zscore"
      )),
      silent = TRUE
    )
  })
}

# saving the plain SIAMCAT model 
saveRDS(sc1, opt$model_rds_out)

# MWMOTE PIPELINE: train, predict, evaluate 
msg("Running MWMOTE pipeline with method: %s", opt$model)
t0_m <- proc.time()[3]
sc_mwm <- siamcat_mwmote_with_custom_ml(
  sc.obj          = sc0,
  opt             = opt,
  minority.class  = "case",
  majority.class  = "control",
  kMinority       = 5,
  seed            = opt$seed
)
sc_mwm <- evaluate.predictions(sc_mwm)
t1_m <- proc.time()[3]

pm_m <- safe_pred_matrix(sc_mwm)
if (is.null(pm_m) || nrow(pm_m) == 0) {
  stop("Prediction matrix is empty in MWMOTE pipeline.")
}
pred_mean_m <- rowMeans(pm_m, na.rm = TRUE)
lbl_m       <- label(sc_mwm)[["label"]]
keep_m <- is.finite(pred_mean_m) & is.finite(lbl_m)
pred_mean_m <- pred_mean_m[keep_m]
lbl_m       <- lbl_m[keep_m]
if (length(pred_mean_m) == 0) stop("All predictions or labels are NA after filtering (MWMOTE).")

yhat_m <- ifelse(pred_mean_m > thr, 1L, 0L)
y_m    <- ifelse(lbl_m == 1L, 1L, 0L)
TP_m <- sum(yhat_m==1L & y_m==1L); TN_m <- sum(yhat_m==0L & y_m==0L)
FP_m <- sum(yhat_m==1L & y_m==0L); FN_m <- sum(yhat_m==0L & y_m==1L)

met_m <- calculate_metrics(TP_m, TN_m, FP_m, FN_m)
met_m$norm   <- opt$norm
met_m$model  <- paste0(opt$model, "_MWMOTE")
met_m$cutoff <- opt$cutoff
met_m$model_id <- sprintf("%s_%s_%s_MWMOTE", opt$norm, opt$model, opt$cutoff)
met_m$train_time_sec <- t1_m - t0_m
met_m$variant <- "mwmote"

# Evaluation plot for MWMOTE predictions only
try(suppressWarnings(model.evaluation.plot(sc_mwm,fn.plot = sub("\\.pdf$", "_mwmote.pdf", opt$eval_pdf))), silent=TRUE)

# save MWMOTE SIAMCAT object if helpful
mwm_rds <- sub("\\.rds$", "_mwmote.rds", opt$model_rds_out)
saveRDS(sc_mwm, mwm_rds)

# write perf rows (plain + mwmote) 
perf_out <- rbind(
  as.data.table(met_plain),
  as.data.table(met_m)
)

setcolorder(perf_out, c("norm","model","cutoff","model_id","variant",
                        "tp","tn","fp","fn","acc","sens","spec","precision",
                        "f1","kappa","mcc","train_time_sec"))

fwrite(perf_out, opt$perf_csv, append = file.exists(opt$perf_csv))

# AUROC curves for both pipelines 
has_pos  <- any(y==1L); has_neg  <- any(y==0L)
nonconst <- (length(unique(pred_mean)) > 1L)
if (has_pos && has_neg && nonconst) {
  roc.obj <- tryCatch(roc(y, pred_mean, quiet=TRUE), error=function(e) NULL)
  if (!is.null(roc.obj)) {
    aucroc_plain <- data.frame(
      FPR           = 1 - roc.obj$specificities,
      TPR           = roc.obj$sensitivities,
      model         = opt$model,
      filter        = opt$cutoff,
      normalization = opt$norm,
      variant       = "siamcat",
      stringsAsFactors = FALSE
    )
  } else {
    aucroc_plain <- NULL
    msg("pROC::roc failed; skipping plain AUROC curve write.")
  }
} else {
  aucroc_plain <- NULL
  msg("Skipping plain AUROC curve (classes or scores degenerate).")
}

# mwmote
has_pos_m  <- any(y_m==1L); has_neg_m  <- any(y_m==0L)
nonconst_m <- (length(unique(pred_mean_m)) > 1L)
if (has_pos_m && has_neg_m && nonconst_m) {
  roc.obj.m <- tryCatch(roc(y_m, pred_mean_m, quiet=TRUE), error=function(e) NULL)
  if (!is.null(roc.obj.m)) {
    aucroc_m <- data.frame(
      FPR           = 1 - roc.obj.m$specificities,
      TPR           = roc.obj.m$sensitivities,
      model         = paste0(opt$model, "_MWMOTE"),
      filter        = opt$cutoff,
      normalization = opt$norm,
      variant       = "mwmote",
      stringsAsFactors = FALSE
    )
  } else {
    aucroc_m <- NULL
    msg("pROC::roc failed; skipping MWMOTE AUROC curve write.")
  }
} else {
  aucroc_m <- NULL
  msg("Skipping MWMOTE AUROC curve (classes or scores degenerate).")
}

auc_combined <- rbind(aucroc_plain, aucroc_m)
if (!is.null(auc_combined) && nrow(auc_combined) > 0) {
  fwrite(auc_combined, opt$auroc_csv, append = file.exists(opt$auroc_csv))
}

msg("Done: %s (plain + MWMOTE)", met_plain$model_id)


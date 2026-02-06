#!/usr/bin/env Rscript

if(!requireNamespace("LiblineaR",quietly=TRUE))
   install.packages("LiblineaR", repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(LiblineaR)
  library(SIAMCAT)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# ---------------- CLI ----------------
opt_list <- list(
  make_option("--siamcat_rds", type="character",
              help="Path to SIAMCAT RDS (contains @norm_feat[['norm.feat']] and labels)"),
  make_option("--model_id",    type="character", default="model",
              help="Model ID used in output filenames"),
  make_option("--method",      type="character", default="ridge_ll",
              help="ridge_ll | lasso_ll | enet | ridge | lasso | randomForest | xgboost"),
  make_option("--mtry",        type="integer",   default=18,
              help="mtry for randomForest (ranger)"),
  make_option("--num_trees",   type="integer",   default=1000,
              help="num.trees for randomForest (ranger)"),
  make_option("--alpha",       type="double",    default=0.5,
              help="alpha for enet (ignored for pure ridge/lasso)"),
  make_option("--threshold",   type="double",    default=0.5,
              help="probability threshold for metrics"),
  make_option("--sample_n",    type="integer",   default=100,
              help="balanced SHAP background size (total rows)"),
  make_option("--nsim",        type="integer",   default=50,
              help="fastshap nsim"),
  make_option("--outdir",      type="character", default=".",
              help="output directory")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$siamcat_rds))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ----------- helpers -----------
canon <- function(s) tolower(gsub("[^a-z0-9]+","", s))
aliases <- list(
  randomforest = "randomForest",
  rf           = "randomForest",
  random_forest= "randomForest",

  ridgell      = "ridge_ll",
  ridgeii      = "ridge_ll",
  ridge_ll     = "ridge_ll",
  ridge_ii     = "ridge_ll",

  lassoll      = "lasso_ll",
  lasso        = "lasso_ll",
  lasso_ll     = "lasso_ll",

  enetll       = "enet",
  enet         = "enet",
  elasticnet   = "enet",
  elasticnetll = "enet",

  xgb          = "xgboost",
  xgboost      = "xgboost"
)
opt$method <- aliases[[canon(opt$method)]] %||% opt$method

pretty_taxon <- function(x) {
  x <- gsub("^(k|p|c|o|f|g|s)__", "", x)
  x <- gsub("\\|", "; ", x)
  x <- gsub("_+", "_", x)
  trimws(x)
}

eval_at <- function(p, y, thr){
  ph <- as.integer(p >= thr)
  TP <- sum(ph==1 & y==1); TN <- sum(ph==0 & y==0)
  FP <- sum(ph==1 & y==0); FN <- sum(ph==0 & y==1)
  ACC  <- (TP+TN)/length(y)
  SENS <- if ((TP+FN)==0) NA_real_ else TP/(TP+FN)
  SPEC <- if ((TN+FP)==0) NA_real_ else TN/(TN+FP)
  data.frame(ACC=ACC,SENS=SENS,SPEC=SPEC,TP=TP,TN=TN,FP=FP,FN=FN)
}

#model training/prediction
fit_predict <- function(method, opt, X, y) {
  # X: matrix/data.frame (n x p), y: 0/1 integer

  if (method == "randomForest") {
    if (!requireNamespace("ranger", quietly=TRUE)) stop("Please install 'ranger'")
    df <- data.frame(label = factor(y), X, check.names = FALSE)
    fit <- ranger::ranger(
      dependent.variable.name = "label",
      data        = df,
      probability = TRUE,
      mtry        = opt$mtry,
      num.trees   = opt$num_trees,
      importance  = "permutation",
      classification = TRUE
    )
    pred_fun <- function(object, newdata) {
      pr <- predict(object, data = data.frame(newdata, check.names = FALSE))$predictions
      as.numeric(pr[, 2])
    }
    return(list(fit=fit, pred_fun=pred_fun))
  }

  if (method %in% c("ridge","lasso","enet")) {
    if (!requireNamespace("glmnet", quietly=TRUE)) stop("Please install 'glmnet'")
    alpha <- if (method == "lasso") 1.0 else if (method == "ridge") 0.0 else opt$alpha
    fit <- glmnet::cv.glmnet(
      x      = as.matrix(X),
      y      = y,
      family = "binomial",
      alpha  = alpha
    )
    pred_fun <- function(object, newdata) {
      as.numeric(predict(object, newx = as.matrix(newdata), s = "lambda.min", type = "response"))
    }
    return(list(fit=fit, pred_fun=pred_fun))
  }

  if (method %in% c("ridge_ll","lasso_ll")) {
    # LiblineaR needs 0/1 integer targets
    type <- if (method == "lasso_ll") 6 else 0
    fit <- LiblineaR::LiblineaR(
      data   = as.matrix(X),
      target = as.integer(y),
      type   = type,
      cost   = 1,
      bias   = 1,
      cross  = 0,
      verbose= FALSE
    )
    pred_fun <- function(object, newdata) {
      pr <- predict(object, as.matrix(newdata), proba = TRUE)
      as.numeric(pr$probabilities[, 2])
    }
    return(list(fit=fit, pred_fun=pred_fun))
  }

  if (method == "xgboost") {
    if (!requireNamespace("xgboost", quietly = TRUE)) stop("Please install 'xgboost'")

    dtrain <- xgboost::xgb.DMatrix(
      data  = as.matrix(X),
      label = as.numeric(y)
    )

    fit <- xgboost::xgboost(
      data      = dtrain,
      objective = "binary:logistic",
      nrounds   = 200,
      verbose   = 0
    )

    pred_fun <- function(object, newdata) {
      # safest: wrap newdata as DMatrix
      dm <- xgboost::xgb.DMatrix(data = as.matrix(newdata))
      as.numeric(predict(object, newdata = dm))
    }

    return(list(fit=fit, pred_fun=pred_fun))
  }

  stop("Unknown --method: ", method)
}

# loading SIAMCAT object
message("[SHAP] loading: ", opt$siamcat_rds)
ssc <- readRDS(opt$siamcat_rds)

# infer variant from filename (for metrics labelling)
variant <- if (grepl("(?i)mwmote", opt$siamcat_rds) ||
              grepl("(?i)mwmote", opt$model_id)) {
  "mwmote"
} else {
  "siamcat"
}
message(sprintf("[SHAP] variant guess: %s", variant))

NF <- tryCatch(get.norm_feat.matrix(ssc), error = function(e) NULL)
if (is.null(NF)) {
  NF <- ssc@norm_feat[["norm.feat"]]
}
stopifnot(is.matrix(NF) || is.data.frame(NF))

feat_names <- rownames(NF)
if (is.null(feat_names)) stop("norm.feat has no rownames; cannot label features.")
feat_names <- trimws(as.character(feat_names))
feat_names[!nzchar(feat_names)] <- paste0("feat_", seq_along(feat_names))
feat_names <- make.unique(feat_names)
feat_pretty <- vapply(feat_names, pretty_taxon, "", USE.NAMES = FALSE)

X <- t(as.matrix(NF))
mode(X) <- "numeric"
X[is.na(X)] <- 0
colnames(X) <- feat_pretty

# Labels (0/1). SIAMCAT convention: 1=case, -1=control
y_raw <- ssc@label$label
y <- as.integer(ifelse(y_raw == -1, 0, 1L))
stopifnot(length(y) == nrow(X))
message(sprintf("[SHAP] data: n=%d p=%d | class0=%d, class1=%d",
                nrow(X), ncol(X), sum(y==0), sum(y==1)))
if (length(unique(y)) < 2) stop("Only one class present; cannot proceed.")

# train model
mp <- fit_predict(opt$method, opt, X, y)
fit <- mp$fit
pred_fun <- mp$pred_fun

# metrics at threshold
p_all <- pred_fun(fit, X)
met <- eval_at(p_all, y, opt$threshold)
met$set      <- "all"
met$method   <- opt$method
met$model_id <- opt$model_id
met$variant  <- variant

met_file <- file.path(opt$outdir, sprintf("metrics_%s_%s.csv", opt$model_id, opt$method))
data.table::fwrite(met, met_file)
message("[SHAP] metrics -> ", met_file)

#  balanced SHAP background 
if (!requireNamespace("fastshap", quietly=TRUE)) stop("Please install 'fastshap'")
if (!requireNamespace("shapviz",  quietly=TRUE)) stop("Please install 'shapviz'")

set.seed(1L)
pos <- which(y==1); neg <- which(y==0)
if (length(pos) == 0 || length(neg) == 0) {
  take <- sample(seq_len(nrow(X)), min(opt$sample_n, nrow(X)))
} else {
  n_each <- max(1L, floor(min(opt$sample_n, nrow(X))/2))
  take <- sort(unique(c(
    sample(pos, min(length(pos), n_each)),
    sample(neg, min(length(neg), n_each))
  )))
}
Xbg <- as.data.frame(X[take, , drop = FALSE])
colnames(Xbg) <- colnames(X)

message(sprintf("[SHAP] background: n=%d (pos=%d,neg=%d) nsim=%d",
                nrow(Xbg), sum(y[take]==1), sum(y[take]==0), opt$nsim))

# SHAP and plots
sv  <- fastshap::explain(
  object       = fit,
  X            = Xbg,
  pred_wrapper = function(object, newdata) pred_fun(object, newdata),
  nsim         = opt$nsim
)
shp <- shapviz::shapviz(sv, X = as.matrix(Xbg))

bar_pdf  <- file.path(opt$outdir,
                      sprintf("%s_%s_%s_SHAP_BAR.pdf",       opt$model_id, opt$method, variant))
bee_pdf  <- file.path(opt$outdir,
                      sprintf("%s_%s_%s_SHAP_Beeswarm.pdf",  opt$model_id, opt$method, variant))
both_pdf <- file.path(opt$outdir,
                      sprintf("%s_%s_%s_SHAP_BAR_BEESWARM.pdf", opt$model_id, opt$method, variant))

grDevices::pdf(bar_pdf, 6, 9)
print(shapviz::sv_importance(shp, max_display = 50))
grDevices::dev.off()

grDevices::pdf(bee_pdf, 6, 9)
print(shapviz::sv_importance(shp, kind = "beeswarm", max_display = 50))
grDevices::dev.off()

grDevices::pdf(both_pdf, 6, 9)
print(shapviz::sv_importance(shp, kind = "both", alpha = 0.2, max_display = 50))
grDevices::dev.off()

message("Done. Wrote: ", paste(c(bar_pdf, bee_pdf, both_pdf), collapse=" | "))


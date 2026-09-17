nextflow.enable.dsl=2

/*******************************
 * Local “read-only” defaults
 * (Populating params directly avoids DSL2 statement conflicts)
 *******************************/
params.features           = params.features           ?: '../Data/FeaturesDiseaseStatusAsiaLung.csv'
params.meta               = params.meta               ?: '../Data/MetaDiseaseStatusAsiaLung.csv'
params.covariates         = params.covariates         ?: 'diversity_shannon,Incidence,Density,Disease_Status'
params.outdir             = params.outdir             ?: '../NextflowResults/results_DiseaseStatusAsia'

params.norms              = params.norms              ?: 'log.std,rank.unit,log.unit'
params.cutoffs            = params.cutoffs            ?: '0.005,0.0001,0.0005,0.01'
params.models             = params.models             ?: 'ridge_ll'
params.rf_consens_thresh  = params.rf_consens_thresh  ?: 0.01

params.sel_metric         = params.sel_metric         ?: 'mcc'
params.balance_lambdas    = params.balance_lambdas    ?: [0.3, 0.4, 0.5]
params.label_column       = params.label_column       ?: 'Disease_Status'
params.case_label         = params.case_label         ?: 'TB_case'
params.taxa_col           = params.taxa_col           ?: 'Genus'
params.meta_id_col        = params.meta_id_col        ?: 'SampleID'

params.shap_nsim          = params.shap_nsim          ?: 100
params.shap_sample        = params.shap_sample        ?: 200
params.shap_mtry          = params.shap_mtry          ?: 'auto'
params.shap_trees         = params.shap_trees         ?: 500
params.shap_alpha         = params.shap_alpha         ?: 0.5
params.shap_thresh        = params.shap_thresh        ?: 0.5

params.val_features       = params.val_features       ?: '../Data/Validation/Validation_Features.csv'
params.val_meta           = params.val_meta           ?: '../Data/Validation/Validation_Meta.csv'
params.val_feat_id_col    = params.val_feat_id_col    ?: params.taxa_col

params.shap_model_ids     = params.shap_model_ids     ?: ''
params.fold_metrics_threshold = params.fold_metrics_threshold ?: 0.5

/*******************************
 * Utils
 *******************************/
def listify(v) {
  if (v == null) return []
  if (v instanceof Collection) return v.collect{ it.toString().trim() }.findAll{ it }
  return v.toString().split(',').collect{ it.trim() }.findAll{ it }
}

/*******************************
 * Processes
 *******************************/
process PREPARE_DATA {
  publishDir "${params.outdir}/base", mode: 'copy'
  input:
    path features
    path meta
    val  seed
  output:
    path "sc_base.rds"
  script:
  """
  Rscript ${projectDir}/bin/prepare_data.R \\
    --features '${features}' \\
    --meta '${meta}' \\
    --taxa_col '${params.taxa_col}' \\
    --meta_id_col '${params.meta_id_col}' \\
    --covars '${params.covariates}' \\
    --label_column '${params.label_column}' \\
    --case_label '${params.case_label}' \\
    --out sc_base.rds
  """
}

process FLT_SPLIT {
  tag "${norm}|${cutoff}"
  publishDir {"${params.outdir}/${norm}"}, mode:'copy'
  input:
    tuple val(norm), val(cutoff), path(base_rds)
  output:
    tuple val(norm),
          val(cutoff),
          path("split_${norm}_${cutoff}.rds"),
          path("association_${cutoff}.pdf"),
          path("confounder_${cutoff}.pdf")
  script:
  """
  Rscript ${projectDir}/bin/flt_split.R \\
    --base ${base_rds} \\
    --norm ${norm} \\
    --cutoff ${cutoff} \\
    --out split_${norm}_${cutoff}.rds \\
    --assoc_pdf association_${cutoff}.pdf \\
    --conf_pdf  confounder_${cutoff}.pdf
  """
}

process TRAIN_EVAL {
  tag "${norm}|${cutoff}|${model}"
  publishDir {"${params.outdir}/${norm}"}, mode:'copy'
  input:
    tuple val(norm), val(cutoff), path(split_rds), val(model)
  output:
    tuple val(norm),
          val(model),
          val(cutoff),
          val("${norm}_${model}_${cutoff}"),
          path("auroc_${norm}_${model}_${cutoff}.csv"),
          path("perf_${norm}_${model}_${cutoff}.csv"),
          path("model_${norm}_${model}_${cutoff}.rds"),
          path("model_${norm}_${model}_${cutoff}_mwmote.rds", optional: true),
          path("evaluation_${norm}_${model}_${cutoff}.pdf", optional: true),
          path("interpretation_${norm}_${model}_${cutoff}.pdf")
  script:
  """
  MODEL_ID=${norm}_${model}_${cutoff}

  Rscript ${projectDir}/bin/train_eval.R \\
    --split ${split_rds} \\
    --model ${model} \\
    --norm ${norm} \\
    --cutoff ${cutoff} \\
    --eval_pdf    evaluation_${norm}_${model}_${cutoff}.pdf \\
    --interp_pdf  interpretation_${norm}_${model}_${cutoff}.pdf \\
    --rf_consens_thresh ${params.rf_consens_thresh} \\
    --auroc_csv   auroc_${norm}_${model}_${cutoff}.csv \\
    --perf_csv    perf_${norm}_${model}_${cutoff}.csv \\
    --model_rds_out model_${norm}_${model}_${cutoff}.rds
  """
}

process FINALIZE_RESULTS {
  publishDir "${params.outdir}", mode:'copy'
  input:
    path all_csvs
  output:
    path "AUROC_plots"
    path "model_performance_merged.csv"
    path "auroc_curve_merged.csv"
  script:
  """
  mkdir -p AUROC_plots

  Rscript ${projectDir}/bin/plot_auroc.R \\
    --indir . \\
    --pattern 'auroc_.*\\.csv' \\
    --out AUROC_plots

  Rscript ${projectDir}/bin/merge_csv.R \\
    --indir . \\
    --perf_out model_performance_merged.csv \\
    --auroc_out auroc_curve_merged.csv
  """
}

process WRITE_MODEL_MAP {
  tag "model-map"
  input:
    tuple val(norm), val(model), val(cutoff), val(model_id),
          path(auroc_csv), path(perf_csv),
          path(model_rds), path(mwmote_rds)
    val outdir_abs
  output:
    path "model_map_${model_id}.tsv"

  script:
  """
  DEST_DIR="${outdir_abs}/${norm}"
  mkdir -p "\${DEST_DIR}"

  PUB_RDS="\${DEST_DIR}/model_${model_id}.rds"
  PUB_RDS_MWM="\${DEST_DIR}/model_${model_id}_mwmote.rds"

  echo -e "model_id\\tmodel_rds" > model_map_${model_id}.tsv
  echo -e "${model_id}\\t\${PUB_RDS}" >> model_map_${model_id}.tsv

  if [ -s "${mwmote_rds}" ]; then
    echo -e "${model_id}_MWMOTE\\t\${PUB_RDS_MWM}" >> model_map_${model_id}.tsv
  fi
  """
}

process FOLD_METRICS {
  tag "fold-metrics"
  publishDir "${params.outdir}", mode:'copy'
  input:
    path model_files
  output:
    path "model_performance_foldwise.csv"
  script:
  """
  set -euo pipefail
  Rscript "${projectDir}/bin/fold_metrics.R" \\
    --indir . \\
    --out model_performance_foldwise.csv \\
    --threshold "${params.fold_metrics_threshold}"
  """
}


process SELECT_TOP3_BY_METRIC {
  tag "select-top3-${params.sel_metric}"
  publishDir "${params.outdir}/top3", mode:'copy'
  input:
    path merged_perf
    path merged_auroc
    path model_maps
    path model_files
  output:
    path "top3.tsv"
  script:
  def lambdas = params.balance_lambdas ?: [0.3, 0.4, 0.5]
  def lambdaStr = lambdas.join(',')
  """
  set -euo pipefail
  cat > select_top3.R <<'RS'
  suppressPackageStartupMessages({ library(data.table) })
  LAMBDAS <- c(${lambdaStr})
  PERF <- fread('model_performance_merged.csv')
  # Filter variant if column exists
  if ("variant" %in% names(PERF)) {
    PERF <- PERF[variant == "siamcat"]
    if (nrow(PERF) == 0L) stop("No rows with variant == 'siamcat' found in model_performance_merged.csv.")
  }
  # Ensure model_id column exists
  if (!'model_id' %in% names(PERF)) {
    if (all(c('norm','model','cutoff') %in% names(PERF))) {
      PERF[, model_id := paste(norm, model, cutoff, sep = '_')]
    } else if ('MODEL_ID' %in% names(PERF)) {
      setnames(PERF, 'MODEL_ID', 'model_id')
    } else stop('model_performance_merged.csv needs model_id or norm/model/cutoff.')
  }
  # Clean string "NA" text and convert target metric columns to numeric
  metric_cols <- c("mcc", "sens", "spec", "auroc")
  for (col in metric_cols) {
    if (col %in% names(PERF)) {
      PERF[get(col) == "NA", (col) := NA]
      PERF[, (col) := as.numeric(as.character(get(col)))]
    }
  }
  # Remove rows where any of the primary sorting metrics are NA
  PERF <- PERF[!is.na(mcc) & !is.na(sens) & !is.na(spec)]
  # Merge with model map files
  map_files <- list.files('.', pattern = glob2rx('model_map_*.tsv'), full.names = TRUE)
  if (!length(map_files)) stop('No model_map_*.tsv files found.')
  MAP <- rbindlist(lapply(map_files, function(f) fread(f, colClasses='character')),
                   fill=TRUE, use.names=TRUE)
  MAP <- unique(MAP[, .(model_id, model_rds)])
  M <- merge(PERF, MAP, by='model_id', all.x=TRUE)
  if (!("auroc" %in% names(M))) M[, auroc := NA_real_]
  M[, gap := abs(sens - spec)]

  # --- Degeneracy guard -----------------------------------------------
  # Flag models whose fitted coefficients are all (near) zero outside the
  # intercept -- this is the failure mode behind "the standard deviation
  # is zero" / constant-prediction warnings seen at validation time. Only
  # applies to coefficient-based models (glmnet ridge/lasso); other model
  # types (e.g. randomForest) are left unflagged (NA) since this check
  # doesn't apply to them.
  is_degenerate <- function(model_id, model_rds, eps = 1e-8) {
    rds_file <- basename(model_rds)
    if (!file.exists(rds_file)) {
      found <- list.files('.', pattern = paste0('^model_', model_id, '\\\\.rds\$'), full.names = TRUE)
      if (!length(found)) return(NA)
      rds_file <- found[1]
    }
    obj <- tryCatch(readRDS(rds_file), error = function(e) NULL)
    if (is.null(obj)) return(NA)
    fit <- obj
    if (is.list(obj) && !is.null(obj\$model)) fit <- obj\$model
    if (!inherits(fit, c("cv.glmnet", "glmnet"))) return(NA)
    co <- tryCatch(as.matrix(coef(fit)), error = function(e) NULL)
    if (is.null(co)) return(NA)
    nm <- rownames(co)
    if (!is.null(nm) && "(Intercept)" %in% nm) co <- co[nm != "(Intercept)", , drop = FALSE]
    isTRUE(all(abs(co) < eps))
  }
  M[, degenerate_coefs := mapply(is_degenerate, model_id, model_rds)]
  n_flagged <- sum(M\$degenerate_coefs, na.rm = TRUE)
  if (n_flagged > 0L) {
    cat(sprintf("[select_top3] Excluding %d model(s) with degenerate (all-zero) coefficients: %s\n",
                n_flagged, paste(M[degenerate_coefs == TRUE, model_id], collapse=", ")))
  }
  M <- M[is.na(degenerate_coefs) | degenerate_coefs == FALSE]
  # ----------------------------------------------------------------------

  # Composite score: MCC penalized by sensitivity/specificity imbalance.
  # Rank across multiple lambda thresholds so the selection isn't dependent
  # on a single choice of penalty weight.
  TOP3 <- rbindlist(lapply(LAMBDAS, function(lam) {
    Mi <- copy(M)
    Mi[, score := mcc - lam * gap]
    setorderv(Mi, cols = "score", order = -1L, na.last = TRUE)
    Mi[1:min(3L, .N), .(lambda = lam,
                        model_id,
                        metric_used = "mcc_minus_lambda_times_sens_spec_gap",
                        metric_value = mcc,
                        sens,
                        spec,
                        gap,
                        score,
                        auroc,
                        model_rds)]
  }))
  fwrite(TOP3, "top3.tsv", sep="\t")
RS
Rscript select_top3.R
  """
}

process SELECT_TOP2_BY_METRIC {
  tag "select-top2-${params.sel_metric}"
  publishDir "${params.outdir}/top2", mode:'copy'
  errorStrategy { task.exitStatus == 1 ? 'ignore' : 'terminate' }
  input:
    path merged_perf
    path merged_auroc
    path model_maps
    path model_files          // NEW: rds files must be staged to check coefficients
  output:
    path "top2.tsv"
  script:
  def lambdas = params.balance_lambdas ?: [0.3, 0.4, 0.5]
  def lambdaStr = lambdas.join(',')
  """
  set -euo pipefail
  cat > select_top2.R <<'RS'
  suppressPackageStartupMessages({ library(data.table) })
  LAMBDAS <- c(${lambdaStr})
  PERF <- fread('model_performance_merged.csv')
  if (!'model_id' %in% names(PERF)) {
    if (all(c('norm','model','cutoff') %in% names(PERF))) {
      PERF[, model_id := paste(norm, model, cutoff, sep = '_')]
    } else if ('MODEL_ID' %in% names(PERF)) {
      setnames(PERF, 'MODEL_ID', 'model_id')
    } else stop('model_performance_merged.csv needs model_id or norm/model/cutoff.')
  }
  metric_cols <- c("mcc", "sens", "spec", "auroc")
  for (col in metric_cols) {
    if (col %in% names(PERF)) {
      PERF[get(col) == "NA", (col) := NA]
      PERF[, (col) := as.numeric(as.character(get(col)))]
    }
  }
  PERF <- PERF[!is.na(mcc) & !is.na(sens) & !is.na(spec)]
  map_files <- list.files('.', pattern = glob2rx('model_map_*.tsv'), full.names = TRUE)
  if (!length(map_files)) stop('No model_map_*.tsv files found.')
  MAP <- rbindlist(lapply(map_files, fread, colClasses='character'), fill=TRUE, use.names=TRUE)
  MAP <- unique(MAP[, .(model_id, model_rds)])
  M <- merge(PERF, MAP, by='model_id', all.x=TRUE)
  if (!("auroc" %in% names(M))) M[, auroc := NA_real_]
  M[, gap := abs(sens - spec)]

  # --- Degeneracy guard (same logic as SELECT_TOP3_BY_METRIC) -----------
  is_degenerate <- function(model_id, model_rds, eps = 1e-8) {
    rds_file <- basename(model_rds)
    if (!file.exists(rds_file)) {
      found <- list.files('.', pattern = paste0('^model_', model_id, '\\\\.rds\$'), full.names = TRUE)
      if (!length(found)) return(NA)
      rds_file <- found[1]
    }
    obj <- tryCatch(readRDS(rds_file), error = function(e) NULL)
    if (is.null(obj)) return(NA)
    fit <- obj
    if (is.list(obj) && !is.null(obj\$model)) fit <- obj\$model
    if (!inherits(fit, c("cv.glmnet", "glmnet"))) return(NA)
    co <- tryCatch(as.matrix(coef(fit)), error = function(e) NULL)
    if (is.null(co)) return(NA)
    nm <- rownames(co)
    if (!is.null(nm) && "(Intercept)" %in% nm) co <- co[nm != "(Intercept)", , drop = FALSE]
    isTRUE(all(abs(co) < eps))
  }
  M[, degenerate_coefs := mapply(is_degenerate, model_id, model_rds)]
  n_flagged <- sum(M\$degenerate_coefs, na.rm = TRUE)
  if (n_flagged > 0L) {
    cat(sprintf("[select_top2] Excluding %d model(s) with degenerate (all-zero) coefficients: %s\n",
                n_flagged, paste(M[degenerate_coefs == TRUE, model_id], collapse=", ")))
  }
  M <- M[is.na(degenerate_coefs) | degenerate_coefs == FALSE]
  # ------------------------------------------------------------------------

  select_top2 <- function(dt, lam) {
    if (nrow(dt) == 0L) return(dt[0])
    di <- copy(dt)
    di[, score := mcc - lam * gap]
    setorderv(di, cols = "score", order = -1L, na.last = TRUE)
    di[1:min(2L, .N)]
  }

  if ("variant" %in% names(M)) {
    M_s <- M[variant == "siamcat"]
    M_m <- M[variant == "mwmote"]
    TOP2 <- rbindlist(lapply(LAMBDAS, function(lam) {
      rbindlist(list(
        select_top2(M_s, lam)[, .(lambda = lam, model_id, variant, mcc, sens, spec, gap, score, auroc, model_rds)],
        select_top2(M_m, lam)[, .(lambda = lam, model_id, variant, mcc, sens, spec, gap, score, auroc, model_rds)]
      ), use.names = TRUE, fill = TRUE)
    }))
  } else {
    TOP2 <- rbindlist(lapply(LAMBDAS, function(lam) {
      select_top2(M, lam)[, .(lambda = lam, model_id, mcc, sens, spec, gap, score, auroc, model_rds)]
    }))
  }
  fwrite(TOP2, 'top2.tsv', sep='\t')
RS
Rscript select_top2.R
  """
}


process VALIDATE_TOP3 {
  tag "validate-top3"
  publishDir "${params.outdir}/validation", mode:'copy'
  input:
    path top3_tsv
    path val_features
    path val_meta
    path model_files
  output:
    path "validation_metrics_*.csv"
    path "validation_evaluation_*.pdf"
    path "validation_roc_*.csv"
  script:
  """
  set -euo pipefail
  # Print top3.tsv content to stdout so it appears in Nextflow task logs for debugging
  echo "=== Contents of ${top3_tsv} ==="
  cat "${top3_tsv}"
  echo "==============================="

  # The file now has a 'lambda' column first and 'model_id' elsewhere in the
  # header (columns are no longer fixed at position 1 / last), and the same
  # model_id can appear multiple times -- once per lambda threshold it was
  # selected under. Look up column positions by header name instead of
  # assuming fixed positions, and validate each unique model_id only once.
  HEADER=\$(head -n1 "${top3_tsv}")
  MODEL_ID_COL=\$(echo "\$HEADER" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="model_id") print i}')
  MODEL_RDS_COL=\$(echo "\$HEADER" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="model_rds") print i}')

  if [ -z "\$MODEL_ID_COL" ] || [ -z "\$MODEL_RDS_COL" ]; then
    echo "ERROR: Could not find 'model_id' and/or 'model_rds' columns in ${top3_tsv} header: \$HEADER" >&2
    exit 1
  fi

  # Deduplicate on model_id: a model selected under more than one lambda
  # only needs to be validated once.
  tail -n +2 "${top3_tsv}" | awk -F'\\t' -v idc="\$MODEL_ID_COL" -v rdsc="\$MODEL_RDS_COL" \\
    '{print \$idc "\\t" \$rdsc}' | sort -u -t\$'\\t' -k1,1 | while IFS=\$'\\t' read -r MODEL_ID MODEL_RDS_RAW; do
    [ -n "\$MODEL_ID" ] || continue
    # Extract base filename
    RDS_FILE=\$(basename "\$MODEL_RDS_RAW")
    # Fallback search if the base filename doesn't match directly
    if [ ! -s "\$RDS_FILE" ]; then
      RDS_FILE=\$(find . -maxdepth 2 -name "model_\${MODEL_ID}.rds" | head -n1)
    fi
    if [ -z "\$RDS_FILE" ] || [ ! -s "\$RDS_FILE" ]; then
      echo "ERROR: Could not find staged model file for ID '\$MODEL_ID' (expected 'model_\${MODEL_ID}.rds')" >&2
      exit 1
    fi
    echo "Running validation for model: \${MODEL_ID} using \${RDS_FILE}..."
    Rscript "${projectDir}/bin/validate.R" \\
      --model_rds    "\$RDS_FILE" \\
      --features     "${val_features}" \\
      --meta         "${val_meta}" \\
      --label_col    "${params.label_column}" \\
      --case_label   "${params.case_label}" \\
      --meta_id_col  "${params.meta_id_col}" \\
      --feat_id_col  "${params.val_feat_id_col}" \\
      --outdir       . \\
      --prefix       "\$MODEL_ID" \\
      --threshold    0.5
  done
  """
}


process SHAP_FROM_SIAMCAT_TOP2 {
  tag "shap-top2"
  publishDir "${params.outdir}/shap_best", mode: 'copy', overwrite: true
  time '6h'
  errorStrategy 'terminate'
  input:
    path top2_tsv
    val  outabs
  output:
    path "*.pdf", optional: true
    path "metrics_*.csv", optional: true
    path "SHAP_RUN.txt", optional: true

  script:
  """
  set -euo pipefail

  run_shap() {
    PROJ="${projectDir}"

    HEADER=\$(head -n1 "${top2_tsv}")
    MODEL_ID_COL=\$(echo "\$HEADER" | awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="model_id") print i}')
    if [ -z "\$MODEL_ID_COL" ]; then
      echo "ERROR: Could not find 'model_id' column in ${top2_tsv} header: \$HEADER" >&2
      exit 1
    fi

    tail -n +2 "${top2_tsv}" | awk -F'\\t' -v idc="\$MODEL_ID_COL" '{print \$idc}' | sort -u | while IFS= read -r MODEL_ID; do
      [ -n "\$MODEL_ID" ] || continue

      NORM=\$(printf "%s" "\$MODEL_ID" | awk -F'_' '{print \$1}')
      CUTOFF=\$(printf "%s" "\$MODEL_ID" | awk -F'_' '{print \$NF}')
      MODEL_RAW=\$(printf "%s" "\$MODEL_ID" | awk -F'_' '{m=""; for(i=2;i<NF;i++){ if(m=="") m=\$i; else m=m"_"\$i } print m }')

      LWR=\$(printf "%s" "\$MODEL_RAW" | tr '[:upper:]' '[:lower:]' | tr '-' '_')
      case "\$LWR" in
        rf|randomforest|random_forest) SHAP_METHOD="randomForest" ;;
        ridge_ll|ridgell|ridge-ll)    SHAP_METHOD="ridge_ll" ;;
        lasso_ll|lasso|lasso-ll)      SHAP_METHOD="lasso_ll" ;;
        enet_ll|enet|elasticnet)      SHAP_METHOD="enet_ll" ;;
        xgboost|xgb)                  SHAP_METHOD="xgboost" ;;
        *)                            SHAP_METHOD="\$MODEL_RAW" ;;
      esac

      RDS="${outabs}/\$NORM/model_\$MODEL_ID.rds"
      if [ ! -s "\$RDS" ]; then
        RDS=\$(find "${outabs}" -type f -name "model_\${MODEL_ID}.rds" -print -quit 2>/dev/null || true)
      fi
      if [ -z "\$RDS" ] || [ ! -s "\$RDS" ]; then
        echo "WARNING: Could not find rds for model_id '\$MODEL_ID' -- skipping." >&2
        continue
      fi

      Rscript "\$PROJ/bin/shap_explain.R" \\
        --siamcat_rds "\$RDS" \\
        --model_id    "\$MODEL_ID" \\
        --method      "\$SHAP_METHOD" \\
        --mtry        "${params.shap_mtry}" \\
        --num_trees   "${params.shap_trees}" \\
        --alpha       "${params.shap_alpha}" \\
        --threshold   "${params.shap_thresh}" \\
        --sample_n    "${params.shap_sample}" \\
        --nsim        "${params.shap_nsim}" \\
        --outdir      .
    done
  }

  run_shap 2>&1 | tee -a SHAP_RUN.txt
  """
}


process SHAP_FROM_MODEL_IDS {
  tag "shap-custom"
  publishDir "${params.outdir}/shap_custom", mode: 'copy', overwrite: true
  time '6h'
  errorStrategy 'terminate'
  input:
    path model_ids_file
    val  outabs
  output:
    path "*.pdf", optional: true
    path "metrics_*.csv", optional: true
    path "SHAP_RUN.txt", optional: true

  script:
  """
  set -euo pipefail

  run_shap() {
    PROJ="${projectDir}"

    while IFS= read -r MODEL_ID; do
      [ -n "\$MODEL_ID" ] || continue

      NORM=\$(printf "%s" "\$MODEL_ID" | awk -F'_' '{print \$1}')
      CUTOFF=\$(printf "%s" "\$MODEL_ID" | awk -F'_' '{print \$NF}')
      MODEL_RAW=\$(printf "%s" "\$MODEL_ID" | awk -F'_' '{m=""; for(i=2;i<NF;i++){ if(m=="") m=\$i; else m=m"_"\$i } print m }')

      LWR=\$(printf "%s" "\$MODEL_RAW" | tr '[:upper:]' '[:lower:]' | tr '-' '_')
      case "\$LWR" in
        rf|randomforest|random_forest) SHAP_METHOD="randomForest" ;;
        ridge_ll|ridgell|ridge-ll)    SHAP_METHOD="ridge_ll" ;;
        lasso_ll|lasso|lasso-ll)      SHAP_METHOD="lasso_ll" ;;
        enet_ll|enet|elasticnet)      SHAP_METHOD="enet_ll" ;;
        xgboost|xgb)                  SHAP_METHOD="xgboost" ;;
        *)                            SHAP_METHOD="\$MODEL_RAW" ;;
      esac

      RDS="${outabs}/\$NORM/model_\$MODEL_ID.rds"
      if [ ! -s "\$RDS" ]; then
        RDS=\$(find "${outabs}" -type f -name "model_\${MODEL_ID}.rds" -print -quit 2>/dev/null || true)
      fi
      if [ -z "\$RDS" ] || [ ! -s "\$RDS" ]; then
        echo "WARNING: Could not find rds for requested model_id '\$MODEL_ID' -- skipping." >&2
        continue
      fi
      RDS_BASENAME=\$(basename "\$RDS")
      case "\$RDS_BASENAME" in
        *_mwmote.rds) SHAP_VARIANT="mwmote" ;;
        *)            SHAP_VARIANT="siamcat" ;;
      esac

      echo "Running SHAP for requested model: \$MODEL_ID (variant: \$SHAP_VARIANT) using \$RDS"

      Rscript "\$PROJ/bin/shap_explain.R" \\
        --siamcat_rds "\$RDS" \\
        --model_id    "\$MODEL_ID" \\
        --method      "\$SHAP_METHOD" \\
        --variant     "\$SHAP_VARIANT" \\
        --mtry        "${params.shap_mtry}" \\
        --num_trees   "${params.shap_trees}" \\
        --alpha       "${params.shap_alpha}" \\
        --threshold   "${params.shap_thresh}" \\
        --sample_n    "${params.shap_sample}" \\
        --nsim        "${params.shap_nsim}" \\
        --outdir      .
    done < "${model_ids_file}"
  }

  run_shap 2>&1 | tee -a SHAP_RUN.txt
  """
}






/*******************************
 * Workflow
 *******************************/
workflow {
  // Move logic inside workflow context to prevent mixing declarations
  def OUTDIR_ABS = file(params.outdir).toAbsolutePath().normalize().toString()
  log.info "OUTDIR_ABS = ${OUTDIR_ABS}"

  // base channels
  CH_FEATURES = Channel.fromPath(params.features)
  CH_META     = Channel.fromPath(params.meta)
  CH_SEED     = Channel.value(42)

  // expand params -> lists
  def NORMS_LIST   = listify(params.norms)
  def CUTOFFS_LIST = listify(params.cutoffs)
  def MODELS_LIST  = listify(params.models)

  def OUTABS = new File(params.outdir).getCanonicalPath()
  CH_OUTABS = Channel.value(OUTABS)

  log.info "NORMS   = ${NORMS_LIST}"
  log.info "CUTOFFS = ${CUTOFFS_LIST}"
  log.info "MODELS  = ${MODELS_LIST}"

  assert NORMS_LIST   && NORMS_LIST.size()   > 0 : "No norms supplied."
  assert CUTOFFS_LIST && CUTOFFS_LIST.size() > 0 : "No cutoffs supplied."
  assert MODELS_LIST  && MODELS_LIST.size()  > 0 : "No models supplied."

  CH_MODELS = Channel.from(MODELS_LIST)

  // 1) prepare once
  PREP = PREPARE_DATA(CH_FEATURES, CH_META, CH_SEED)

  // 2) build (norm, cutoff) grid
  def COMBO_LIST = NORMS_LIST.collectMany { n -> CUTOFFS_LIST.collect { c -> tuple(n, c) } }
  CH_COMBO = Channel.fromList(COMBO_LIST)

  // 3) pair base with each combo
  SPLIT_IN = CH_COMBO.combine(PREP).map { it -> tuple(it[0], it[1], it[2]) }

  // 4) split/filter
  SPLITS = FLT_SPLIT(SPLIT_IN)

  // 5) cross with models
  GRID = SPLITS.combine(CH_MODELS).map { z -> tuple(z[0], z[1], z[2], z[5]) }

  // 6) train/eval
  TRAINED = TRAIN_EVAL(GRID)

  // 7) Collect all RDS files created during training
   //ALL_MODEL_RDS = TRAINED.map { t -> t[6] }.collect()
   ALL_MODEL_RDS = TRAINED.map { t -> [t[6], t[7]] }.flatten().collect()
   //   FOLDWISE_CSV = FOLD_METRICS(ALL_MODEL_RDS)

  // 8) gather perf + auroc csvs
  ALL_CSVS = TRAINED.map { t -> [t[4], t[5]] }.flatten().collect()

  // 9) finalize once
  def (CH_AUROC_PLOTS, CH_MODEL_PERF_MERGED, CH_AUROC_MERGED) = FINALIZE_RESULTS(ALL_CSVS)

  // 10) model maps
  MODEL_MAP_SHARDS = TRAINED.map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]) }
  MODEL_MAP_FILES  = WRITE_MODEL_MAP(MODEL_MAP_SHARDS, CH_OUTABS).collect()

   // 11) selections & validation
  TOP3_TSV = SELECT_TOP3_BY_METRIC(CH_MODEL_PERF_MERGED, CH_AUROC_MERGED, MODEL_MAP_FILES, ALL_MODEL_RDS)
  VALIDATE_TOP3(TOP3_TSV, Channel.fromPath(params.val_features), Channel.fromPath(params.val_meta), ALL_MODEL_RDS)
  // 12) SHAP from top 2
  TOP2_TSV = SELECT_TOP2_BY_METRIC(CH_MODEL_PERF_MERGED, CH_AUROC_MERGED, MODEL_MAP_FILES, ALL_MODEL_RDS)
  SHAP_FROM_SIAMCAT_TOP2(TOP2_TSV, CH_OUTABS)


  // 13) SHAP on explicitly requested model_ids
   def SHAP_MODEL_IDS_LIST = listify(params.shap_model_ids)
    if (SHAP_MODEL_IDS_LIST) {
     log.info "SHAP_MODEL_IDS = ${SHAP_MODEL_IDS_LIST}"
    CH_SHAP_IDS_FILE = Channel
      .fromList(SHAP_MODEL_IDS_LIST)
      .collectFile(name: 'requested_model_ids.txt', newLine: true)
    SHAP_FROM_MODEL_IDS(CH_SHAP_IDS_FILE, CH_OUTABS)
    }


}

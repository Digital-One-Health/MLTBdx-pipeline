nextflow.enable.dsl=2

/*******************************
 * Local “read-only” defaults
 * (Define real defaults in nextflow.config; these are just fallbacks)
 *******************************/
def FEATURES   = params.features   ?: '../Data/FeaturesDiseaseStatusAsiaLung.csv'
def META       = params.meta       ?: '../Data/MetaDiseaseStatusAsiaLung.csv'
def COVARIATES = params.covariates  ?: 'diversity_shannon,Incidence,Density,Disease_Status'
def OUTDIR     = params.outdir     ?: '../NextflowResults/results_DiseaseStatusAsia'

def NORMS_RAW  = params.norms      ?: 'log.std,rank.unit,log.unit'
def CUTOFFS_RAW= params.cutoffs    ?: '0.005,0.0001,0.0005,0.01'
def MODELS_RAW = params.models     ?: 'ridge_ll'
def RF_CONSENS_THRESH = (params.rf_consens_thresh   ?: 0.01) as float

def SEL_METRIC = params.sel_metric ?: 'mcc'               // single canonical name
def LABEL_COL  = (params.label_column ?: params.label_col) ?: 'Disease_Status'
def CASE_LABEL = params.case_label   ?: 'TB_case'
def TAXA_COL   = params.taxa_col     ?: 'Genus'
def META_IDCOL = params.meta_id_col  ?: 'SampleID'

def SHAP_NSIM  = (params.shap_nsim   ?: 100)  as int
def SHAP_SAMP  = (params.shap_sample ?: 200)  as int
def SHAP_MTRY  = (params.shap_mtry   ?: 18)   as int
def SHAP_TREES = (params.shap_trees  ?: 1000) as int
def SHAP_ALPHA = (params.shap_alpha  ?: 0.5)  as float
def SHAP_THRESH= (params.shap_thresh ?: 0.5)  as float

def VAL_FEATURES    = params.val_features    ?: '../Data/Validation/Validation_Features.csv'
def VAL_META        = params.val_meta        ?: '../Data/Validation/Validation_Meta.csv'
def VAL_FEAT_IDCOL  = params.val_feat_id_col ?: TAXA_COL

/*******************************
 * Utils
 *******************************/
def listify(v) {
  if (v == null) return []
  if (v instanceof Collection) return v.collect{ it.toString().trim() }.findAll{ it }
  return v.toString().split(',').collect{ it.trim() }.findAll{ it }
}

// Make outdir absolute and canonical
def OUTDIR_ABS = file(OUTDIR).toAbsolutePath().normalize().toString()
log.info "OUTDIR_ABS = ${OUTDIR_ABS}"

/*******************************
 * Processes
 *******************************/
process PREPARE_DATA {
  publishDir "${OUTDIR}/base", mode: 'copy'
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
    --taxa_col '${TAXA_COL}' \\
    --meta_id_col '${META_IDCOL}' \\
    --covars '${COVARIATES}' \\
    --label_column '${LABEL_COL}' \\
    --case_label '${CASE_LABEL}' \\
    --out sc_base.rds
  """
}

process FLT_SPLIT {
  tag "${norm}|${cutoff}"
  publishDir "${OUTDIR}/${norm}", mode:'copy'
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
  publishDir "${OUTDIR}/${norm}", mode:'copy'
  input:
    tuple val(norm), val(cutoff), path(split_rds), val(model)
  output:
    tuple val(norm),
          val(model),
          val(cutoff),
          val("${norm}_${model}_${cutoff}"),                 // model_id (plain SIAMCAT)
          path("auroc_${norm}_${model}_${cutoff}.csv"),
          path("perf_${norm}_${model}_${cutoff}.csv"),
          path("model_${norm}_${model}_${cutoff}.rds"),
          path("model_${norm}_${model}_${cutoff}_mwmote.rds", optional: true),   // ### NEW: MWMOTE RDS
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
    --rf_consens_thresh ${RF_CONSENS_THRESH} \\
    --auroc_csv   auroc_${norm}_${model}_${cutoff}.csv \\
    --perf_csv    perf_${norm}_${model}_${cutoff}.csv \\
    --model_rds_out model_${norm}_${model}_${cutoff}.rds
  """
}

process FINALIZE_RESULTS {
  publishDir "${OUTDIR}", mode:'copy'
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
  publishDir "${OUTDIR_ABS}/model_maps", mode: 'copy'

  input:
    tuple val(norm), val(model), val(cutoff), val(model_id),
          path(auroc_csv), path(perf_csv),
          path(model_rds),                   // plain SIAMCAT
          path(mwmote_rds)   // ### NEW: MWMOTE RDS
  output:
    path "model_map_${model_id}.tsv"

  script:
  """
  DEST_DIR="${OUTDIR_ABS}/${norm}"
  mkdir -p "\${DEST_DIR}"

  PUB_RDS="\${DEST_DIR}/model_${model_id}.rds"
  PUB_RDS_MWM="\${DEST_DIR}/model_${model_id}_mwmote.rds"

  echo -e "model_id\\tmodel_rds" > model_map_${model_id}.tsv
  echo -e "${model_id}\\t\${PUB_RDS}" >> model_map_${model_id}.tsv

  ### NEW: add mapping for the MWMOTE variant if the file exists
  if [ -s "${mwmote_rds}" ]; then
    echo -e "${model_id}_MWMOTE\\t\${PUB_RDS_MWM}" >> model_map_${model_id}.tsv
  fi
  ### END NEW
  """
}


process SELECT_TOP3_BY_METRIC {
  tag "select-top3-${SEL_METRIC}"
  publishDir "${OUTDIR}/top3", mode:'copy'

  input:
    path merged_perf
    path merged_auroc
    path model_maps

  output:
    path "top3.tsv"

  script:
  """
  set -euo pipefail

  cat > select_top3.R <<'RS'
  suppressPackageStartupMessages({ library(data.table) })

  PERF <- fread('model_performance_merged.csv')

  # Only SIAMCAT models are valid for external validation
  if ("variant" %in% names(PERF)) {
    PERF <- PERF[variant == "siamcat"]
    if (nrow(PERF) == 0L) stop("No rows with variant == 'siamcat' found in model_performance_merged.csv.")
  }

  # Ensure model_id exists
  if (!'model_id' %in% names(PERF)) {
    if (all(c('norm','model','cutoff') %in% names(PERF))) {
      PERF[, model_id := paste(norm, model, cutoff, sep = '_')]
    } else if ('MODEL_ID' %in% names(PERF)) {
      setnames(PERF, 'MODEL_ID', 'model_id')
    } else stop('model_performance_merged.csv needs model_id or norm/model/cutoff.')
  }

  # Try to standardize AUROC column if present; otherwise create placeholder
  cand_auc <- c('auroc','AUC','AUROC','auc')
  auc_col <- cand_auc[cand_auc %in% names(PERF)][1]
  if (!is.na(auc_col)) {
    setnames(PERF, auc_col, 'auroc')
  } else {
    PERF[, auroc := NA_real_]
  }

  # column matching (case-insensitive, alnum normalized)
  desired <- tolower("${SEL_METRIC}")
  alnum <- function(x) gsub('[^a-z0-9]+','', tolower(x))
  cn_alnum <- alnum(names(PERF))
  metric_idx <- which(cn_alnum == alnum(desired))[1]
  metric_col <- if (length(metric_idx)) names(PERF)[metric_idx] else NA_character_

  # If desired metric missing, fall back
  if (is.na(metric_col) || !(metric_col %in% names(PERF))) {
    metric_col <- if ("mcc" %in% names(PERF)) "mcc" else if ("accuracy" %in% names(PERF)) "accuracy" else "auroc"
  }
  if (!(metric_col %in% names(PERF))) metric_col <- "model_id"

  # Read model maps written by WRITE_MODEL_MAP
  map_files <- list.files('.', pattern = glob2rx('model_map_*.tsv'), full.names = TRUE)
  if (!length(map_files)) stop('No model_map_*.tsv files found.')

  MAP <- rbindlist(lapply(map_files, function(f) fread(f, colClasses='character')),
                   fill=TRUE, use.names=TRUE)
  MAP <- unique(MAP[, .(model_id, model_rds)])

  # Merge perf + map
  M <- merge(PERF, MAP, by='model_id', all.x=TRUE)

  # Ensure auroc exists on merged too
  if (!("auroc" %in% names(M))) M[, auroc := NA_real_]

  # Sort only by columns that exist
  sort_cols <- unique(c(metric_col, "auroc", "accuracy", "model_id"))
  sort_cols <- sort_cols[sort_cols %in% names(M)]

  ord <- rep(-1L, length(sort_cols))
  if (length(sort_cols) && tail(sort_cols, 1L) == "model_id") ord[length(ord)] <- +1L

  if (length(sort_cols)) setorderv(M, sort_cols, ord, na.last = TRUE)

  # Take top 3
  TOP3 <- M[1:min(3L, .N), .(model_id,
                            metric_used = metric_col,
                            metric_value = get(metric_col),
                            auroc,
                            model_rds)]
  fwrite(TOP3, "top3.tsv", sep="\\t")

  RS

  Rscript select_top3.R
  """
}



process SELECT_TOP2_BY_METRIC {
  tag "select-top2-${SEL_METRIC}"
  publishDir "${OUTDIR}/top2", mode:'copy'

  input:
    path merged_perf
    path merged_auroc
    path model_maps

  output:
    path "top2.tsv"

  errorStrategy {
    task.exitStatus == 1 ? 'ignore' : 'terminate'
  }

  script:
  """
  set -euo pipefail
  cat > select_top2.R <<'RS'
  suppressPackageStartupMessages({ library(data.table) })
  PERF <- fread('model_performance_merged.csv')

  if (!'model_id' %in% names(PERF)) {
    if (all(c('norm','model','cutoff') %in% names(PERF))) {
      PERF[, model_id := paste(norm, model, cutoff, sep = '_')]
    } else if ('MODEL_ID' %in% names(PERF)) {
      setnames(PERF, 'MODEL_ID', 'model_id')
    } else stop('model_performance_merged.csv needs model_id or norm/model/cutoff.')
  }

  cand_auc <- c('auroc','AUC','AUROC','auc')
  auc_in_perf <- cand_auc[cand_auc %in% names(PERF)][1]
  if (!is.na(auc_in_perf)) setnames(PERF, auc_in_perf, 'auroc') else PERF[, auroc := NA_real_]

  desired <- tolower("${SEL_METRIC}")
  alias <- list('f1'='f1','f1_score'='f1','f1score'='f1','mcc'='mcc',
                'accuracy'='accuracy','sens'='sensitivity','recall'='sensitivity','tpr'='sensitivity',
                'spec'='specificity','tnr'='specificity','prec'='precision','ppv'='precision',
                'kappa'='kappa','auroc'='auroc','auc'='auroc')
  key <- if (desired %in% names(alias)) alias[[desired]] else desired
  alnum <- function(x) gsub('[^a-z0-9]+','', tolower(x))
  cn_alnum <- alnum(names(PERF))
  target_idx <- which(cn_alnum == alnum(key))[1]
  metric_col <- if (length(target_idx)) names(PERF)[target_idx] else NA_character_
  if (is.na(metric_col) || !(metric_col %in% names(PERF))) {
    warning(sprintf("Metric '%s' not found; falling back to AUROC.", key))
    metric_col <- 'auroc'
  }

  map_files <- list.files('.', pattern = glob2rx('model_map_*.tsv'), full.names = TRUE)
  if (!length(map_files)) stop('No model_map_*.tsv files found.')
  MAP <- rbindlist(lapply(map_files, fread, colClasses='character'), fill=TRUE, use.names=TRUE)
  MAP <- unique(MAP[, .(model_id, model_rds)])

  M <- merge(PERF, MAP, by='model_id', all.x=TRUE)

  sort_cols <- c(metric_col, 'auroc', 'accuracy', 'model_id')
  sort_cols <- sort_cols[sort_cols %in% names(M)]
  ord <- rep(-1L, length(sort_cols))
  if (length(sort_cols) && tail(sort_cols,1L) == 'model_id') ord[length(ord)] <- +1L
  if (length(sort_cols)) data.table::setorderv(M, sort_cols, ord, na.last=TRUE)

  # top-2 SIAMCAT + top-2 MWMOTE
  if ("variant" %in% names(M)) {
    M_s <- M[variant == "siamcat"]
    M_m <- M[variant == "mwmote"]

    TOP_S <- if (nrow(M_s) > 0) M_s[1:min(2L, .N), .(model_id, auroc, model_rds)] else M[0]
    TOP_M <- if (nrow(M_m) > 0) M_m[1:min(2L, .N), .(model_id, auroc, model_rds)] else M[0]

    TOP2 <- rbind(TOP_S, TOP_M, use.names = TRUE, fill = TRUE)
  } else {
    TOP2 <- M[1:min(2L, .N), .(model_id, auroc, model_rds)]
  }
  ### END NEW

  fwrite(TOP2, 'top2.tsv', sep='\\t')
  RS

  Rscript select_top2.R
  """
}



process VALIDATE_TOP3 {
  tag "validate-top3"
  publishDir "${OUTDIR}/validation", mode:'copy'

  input:
    path top3_tsv
    path val_features
    path val_meta

  output:
    path "validation_metrics_*.csv"
    path "validation_evaluation_*.pdf"
    path "validation_roc_*.csv"

  script:
  """
  set -euo pipefail

  # Use staged filenames (do NOT reference bash variables like \$top3_tsv)
  TOP3_FILE="${top3_tsv}"
  VAL_FEAT="${val_features}"
  VAL_META="${val_meta}"

  echo "[VALIDATE] top3 file: \$TOP3_FILE"
  echo "[VALIDATE] val feat : \$VAL_FEAT"
  echo "[VALIDATE] val meta : \$VAL_META"
  ls -l

  # Read TSV reliably with cut (no IFS=\$'\\t' problems)
  tail -n +2 "\$TOP3_FILE" | while read -r LINE; do
    MODEL_ID=\$(printf "%s" "\$LINE" | cut -f1)
    MODEL_RDS=\$(printf "%s" "\$LINE" | cut -f5)

    echo "[VALIDATE] Model: \$MODEL_ID"
    echo "[VALIDATE] RDS:   \$MODEL_RDS"

    [[ -n "\$MODEL_ID" && -n "\$MODEL_RDS" ]] || { echo "[VALIDATE] Bad row: \$LINE" >&2; continue; }
    [[ -s "\$MODEL_RDS" ]] || { echo "[VALIDATE] Missing RDS: \$MODEL_RDS" >&2; continue; }

    Rscript "${projectDir}/bin/validate.R" \\
      --model_rds   "\$MODEL_RDS" \\
      --features    "\$VAL_FEAT" \\
      --meta        "\$VAL_META" \\
      --label_col   "${LABEL_COL}" \\
      --case_label  "${CASE_LABEL}" \\
      --meta_id_col "${META_IDCOL}" \\
      --feat_id_col "${VAL_FEAT_IDCOL}" \\
      --outdir      . \\
      --prefix      "\$MODEL_ID" \\
      --threshold   0.5
  done
  """
}


process SHAP_FROM_SIAMCAT_TOP2 {
  tag "shap-top2"
  publishDir "${OUTDIR}/shap_best", mode: 'copy', overwrite: true
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
  exec > >(stdbuf -oL -eL tee -a SHAP_RUN.txt) 2>&1

  PROJ="${projectDir}"

  echo "[SHAP] top2: ${top2_tsv}"
  echo "[SHAP] results root(abs): ${outabs}"
  ls -1 "${outabs}" || true

  tail -n +2 "${top2_tsv}" | cut -f1 | while IFS= read -r MODEL_ID; do
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
      *)                             SHAP_METHOD="\$MODEL_RAW" ;;
    esac

    RDS="${outabs}/\$NORM/model_\$MODEL_ID.rds"
    echo "[SHAP] probe: \$RDS"
    if [ ! -s "\$RDS" ]; then
      echo "[SHAP] not found, searching under ${outabs}"
      RDS=\$(find "${outabs}" -type f -name "model_\${MODEL_ID}.rds" -print -quit 2>/dev/null || true)
    fi
    [ -n "\$RDS" ] && [ -s "\$RDS" ] || { echo "[SHAP] missing RDS for \$MODEL_ID"; continue; }

    echo "[SHAP] MODEL_ID=\$MODEL_ID | method=\$SHAP_METHOD"
    echo "[SHAP] RDS=\$RDS"

    Rscript "\$PROJ/bin/shap_explain.R" \\
      --siamcat_rds "\$RDS" \\
      --model_id    "\$MODEL_ID" \\
      --method      "\$SHAP_METHOD" \\
      --mtry        "${SHAP_MTRY}" \\
      --num_trees   "${SHAP_TREES}" \\
      --alpha       "${SHAP_ALPHA}" \\
      --threshold   "${SHAP_THRESH}" \\
      --sample_n    "${SHAP_SAMP}" \\
      --nsim        "${SHAP_NSIM}" \\
      --outdir      .
  done

  echo "[SHAP] done."
  """
}

/*******************************
 * Workflow
 *******************************/
workflow {
  // base channels
  CH_FEATURES = Channel.fromPath(FEATURES)
  CH_META     = Channel.fromPath(META)
  CH_SEED     = Channel.value(42)

  // expand params -> lists
  def NORMS_LIST   = listify(NORMS_RAW)
  def CUTOFFS_LIST = listify(CUTOFFS_RAW)
  def MODELS_LIST  = listify(MODELS_RAW)

  // Make absolute path to your results root once
  def OUTABS = new File(OUTDIR).getCanonicalPath()
  CH_OUTABS = Channel.value(OUTABS)

  log.info "NORMS   = ${NORMS_LIST}"
  log.info "CUTOFFS = ${CUTOFFS_LIST}"
  log.info "MODELS  = ${MODELS_LIST}"

  assert NORMS_LIST   && NORMS_LIST.size()   > 0 : "No norms supplied. Try --norms 'log.std,rank.unit,log.unit'"
  assert CUTOFFS_LIST && CUTOFFS_LIST.size() > 0 : "No cutoffs supplied. Try --cutoffs '0.005,0.0001,0.0005,0.01'"
  assert MODELS_LIST  && MODELS_LIST.size()  > 0 : "No models supplied. Try --models 'ridge_ll'"

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

  // 7) gather perf + auroc csvs
  ALL_CSVS = TRAINED.map { t -> [t[4], t[5]] }.flatten().collect()

  // 8) finalize once
  def (CH_AUROC_PLOTS, CH_MODEL_PERF_MERGED, CH_AUROC_MERGED) = FINALIZE_RESULTS(ALL_CSVS)

  // 9) model maps (now include both SIAMCAT + MWMOTE RDS mapping)
  MODEL_MAP_SHARDS = TRAINED.map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]) }   // ### NEW: pass mwmote_rds
  MODEL_MAP_FILES  = WRITE_MODEL_MAP(MODEL_MAP_SHARDS).collect()

  // 10) selections & validation
  TOP3_TSV = SELECT_TOP3_BY_METRIC(CH_MODEL_PERF_MERGED, CH_AUROC_MERGED, MODEL_MAP_FILES)
  VALIDATE_TOP3(TOP3_TSV, Channel.fromPath(VAL_FEATURES), Channel.fromPath(VAL_META))


  TOP2_TSV = SELECT_TOP2_BY_METRIC(CH_MODEL_PERF_MERGED, CH_AUROC_MERGED, MODEL_MAP_FILES)
  SHAP_FROM_SIAMCAT_TOP2(TOP2_TSV, CH_OUTABS)
}

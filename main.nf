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
params.label_column       = params.label_column       ?: 'Disease_Status'
params.case_label         = params.case_label         ?: 'TB_case'
params.taxa_col           = params.taxa_col           ?: 'Genus'
params.meta_id_col        = params.meta_id_col        ?: 'SampleID'

params.shap_nsim          = params.shap_nsim          ?: 100
params.shap_sample        = params.shap_sample        ?: 200
params.shap_mtry          = params.shap_mtry          ?: 18
params.shap_trees         = params.shap_trees         ?: 1000
params.shap_alpha         = params.shap_alpha         ?: 0.5
params.shap_thresh        = params.shap_thresh        ?: 0.5

params.val_features       = params.val_features       ?: '../Data/Validation/Validation_Features.csv'
params.val_meta           = params.val_meta           ?: '../Data/Validation/Validation_Meta.csv'
params.val_feat_id_col    = params.val_feat_id_col    ?: params.taxa_col

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

process SELECT_TOP3_BY_METRIC {
  tag "select-top3-${params.sel_metric}"
  publishDir "${params.outdir}/top3", mode:'copy'
  
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

  # Apply multi-level sort: MCC (1st), SENS (2nd), SPEC (3rd) — all descending
  setorderv(M, cols = c("mcc", "sens", "spec"), order = c(-1L, -1L, -1L), na.last = TRUE)

  #  Extract Top 3 models
  if (!("auroc" %in% names(M))) M[, auroc := NA_real_]

  TOP3 <- M[1:min(3L, .N), .(model_id,
                             metric_used = "mcc_sens_spec",
                             metric_value = mcc,
                             sens,
                             spec,
                             auroc,
                             model_rds)]

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

  output:
    path "top2.tsv"

  script:
  """
  set -euo pipefail

  cat > select_top2.R <<'RS'
  suppressPackageStartupMessages({ library(data.table) })

  PERF <- fread('model_performance_merged.csv')

  # Ensure model_id exists
  if (!'model_id' %in% names(PERF)) {
    if (all(c('norm','model','cutoff') %in% names(PERF))) {
      PERF[, model_id := paste(norm, model, cutoff, sep = '_')]
    } else if ('MODEL_ID' %in% names(PERF)) {
      setnames(PERF, 'MODEL_ID', 'model_id')
    } else stop('model_performance_merged.csv needs model_id or norm/model/cutoff.')
  }

  # Clean string "NA" text and coerce target metrics to numeric
  metric_cols <- c("mcc", "sens", "spec", "auroc")
  for (col in metric_cols) {
    if (col %in% names(PERF)) {
      PERF[get(col) == "NA", (col) := NA]
      PERF[, (col) := as.numeric(as.character(get(col)))]
    }
  }

  # Filter out rows missing core sorting metrics
  PERF <- PERF[!is.na(mcc) & !is.na(sens) & !is.na(spec)]

  # Merge with model map files
  map_files <- list.files('.', pattern = glob2rx('model_map_*.tsv'), full.names = TRUE)
  if (!length(map_files)) stop('No model_map_*.tsv files found.')

  MAP <- rbindlist(lapply(map_files, fread, colClasses='character'), fill=TRUE, use.names=TRUE)
  MAP <- unique(MAP[, .(model_id, model_rds)])

  M <- merge(PERF, MAP, by='model_id', all.x=TRUE)

  # Apply multi-level sort: MCC (1st), SENS (2nd), SPEC (3rd) — all descending
  setorderv(M, cols = c("mcc", "sens", "spec"), order = c(-1L, -1L, -1L), na.last = TRUE)

  # Select Top 2 per variant (if present) or overall Top 2
  if (!("auroc" %in% names(M))) M[, auroc := NA_real_]

  if ("variant" %in% names(M)) {
    M_s <- M[variant == "siamcat"]
    M_m <- M[variant == "mwmote"]

    TOP_S <- if (nrow(M_s) > 0) M_s[1:min(2L, .N), .(model_id, variant, mcc, sens, spec, auroc, model_rds)] else M[0]
    TOP_M <- if (nrow(M_m) > 0) M_m[1:min(2L, .N), .(model_id, variant, mcc, sens, spec, auroc, model_rds)] else M[0]

    TOP2 <- rbind(TOP_S, TOP_M, use.names = TRUE, fill = TRUE)
  } else {
    TOP2 <- M[1:min(2L, .N), .(model_id, mcc, sens, spec, auroc, model_rds)]
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
  output:
    path "validation_metrics_*.csv"
    path "validation_evaluation_*.pdf"
    path "validation_roc_*.csv"

  script:
  """
  set -euo pipefail

  TOP3_FILE="${top3_tsv}"
  VAL_FEAT="${val_features}"
  VAL_META="${val_meta}"

  tail -n +2 "\$TOP3_FILE" | while read -r LINE; do
    MODEL_ID=\$(printf "%s" "\$LINE" | cut -f1)
    MODEL_RDS=\$(printf "%s" "\$LINE" | cut -f5)

    [[ -n "\$MODEL_ID" && -n "\$MODEL_RDS" ]] || continue
    [[ -s "\$MODEL_RDS" ]] || continue

    Rscript "${projectDir}/bin/validate.R" \\
      --model_rds   "\$MODEL_RDS" \\
      --features    "\$VAL_FEAT" \\
      --meta        "\$VAL_META" \\
      --label_col   "${params.label_column}" \\
      --case_label  "${params.case_label}" \\
      --meta_id_col "${params.meta_id_col}" \\
      --feat_id_col "${params.val_feat_id_col}" \\
      --outdir      . \\
      --prefix      "\$MODEL_ID" \\
      --threshold   0.5
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
  exec > >(stdbuf -oL -eL tee -a SHAP_RUN.txt) 2>&1

  PROJ="${projectDir}"

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
      *)                            SHAP_METHOD="\$MODEL_RAW" ;;
    esac

    RDS="${outabs}/\$NORM/model_\$MODEL_ID.rds"
    if [ ! -s "\$RDS" ]; then
      RDS=\$(find "${outabs}" -type f -name "model_\${MODEL_ID}.rds" -print -quit 2>/dev/null || true)
    fi
    [ -n "\$RDS" ] && [ -s "\$RDS" ] || continue

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

  // 7) gather perf + auroc csvs
  ALL_CSVS = TRAINED.map { t -> [t[4], t[5]] }.flatten().collect()

  // 8) finalize once
  def (CH_AUROC_PLOTS, CH_MODEL_PERF_MERGED, CH_AUROC_MERGED) = FINALIZE_RESULTS(ALL_CSVS)

  // 9) model maps 
  MODEL_MAP_SHARDS = TRAINED.map { t -> tuple(t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]) } 
  MODEL_MAP_FILES  = WRITE_MODEL_MAP(MODEL_MAP_SHARDS, CH_OUTABS).collect()

  // 10) selections & validation
  TOP3_TSV = SELECT_TOP3_BY_METRIC(CH_MODEL_PERF_MERGED, CH_AUROC_MERGED, MODEL_MAP_FILES)
  VALIDATE_TOP3(TOP3_TSV, Channel.fromPath(params.val_features), Channel.fromPath(params.val_meta))

  TOP2_TSV = SELECT_TOP2_BY_METRIC(CH_MODEL_PERF_MERGED, CH_AUROC_MERGED, MODEL_MAP_FILES)
  SHAP_FROM_SIAMCAT_TOP2(TOP2_TSV, CH_OUTABS)
}

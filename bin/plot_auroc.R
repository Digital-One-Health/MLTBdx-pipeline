#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(ggplot2); library(dplyr); library(readr); library(stringr)
})

opt_list <- list(
  make_option("--indir",   type="character", default="."),
  make_option("--pattern", type="character", default="auroc_.*\\.csv"),
  make_option("--out",     type="character", default="AUROC_plots")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# Find files
files <- list.files(opt$indir, pattern = opt$pattern, full.names = TRUE, recursive = TRUE)
if (length(files) == 0) {
  message("No AUROC CSV files found in: ", opt$indir, " with pattern: ", opt$pattern)
  quit(save="no", status=0)
}

# Read, combine and normalization
read_one <- function(f) {
  df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
  # Standardize column names
  names(df) <- gsub("\\s+", "", names(df))
  # Ensure required columns exist
  need <- c("FPR","TPR","model","filter")
  if (!all(need %in% names(df))) {
    stop("File lacks required columns FPR, TPR, model, filter: ", f)
  }
  if (!("normalization" %in% names(df))) {
    # try infer from path
    norm <- str_match(f, "/(log\\.std|rank\\.unit|log\\.unit)/")[,2]
    if (is.na(norm)) norm <- "merged"
    df$normalization <- norm
  }
  df$source_file <- basename(f)
  df
}

au <- bind_rows(lapply(files, read_one))

# Ensuring the right datatypes
au <- au %>%
  mutate(
    model  = as.factor(model),
    filter = as.factor(filter),
    normalization = as.factor(normalization)
  )

# Function nice ROC plot for a given data frame
plot_roc <- function(df, title_suffix = "") {
  df <- df %>% arrange(model, filter, FPR, TPR)

  ggplot(df, aes(x = FPR, y = TPR,
                 group = interaction(model, filter),
                 color = filter)) +
    geom_path(linewidth = 0.8, alpha = 0.9, na.rm = TRUE) +
    facet_wrap(~ model, ncol = 2, scales = "fixed") +   # <- was "free_y"
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +       # square ROC box
    theme_minimal(base_size = 12) +
    labs(title = paste0("AUROC Curves by Filter & Model", title_suffix),
         x = "False Positive Rate", y = "True Positive Rate", color = "Filter")
}

# 1) Combined plot over all data
p_all <- plot_roc(au, " — All")
ggsave(file.path(opt$out, "auroc_ALL.pdf"), plot = p_all, width = 11, height = 7)

# 2) One PDF per normalization level
for (nm in levels(au$normalization)) {
  df <- au %>% filter(normalization == nm)
  if (nrow(df) == 0) next
  p <- plot_roc(df, paste0(" — ", nm))
  ggsave(file.path(opt$out, paste0("auroc_", nm, ".pdf")), plot = p, width = 11, height = 7)
}

message("Done. Wrote PDFs to: ", opt$out)


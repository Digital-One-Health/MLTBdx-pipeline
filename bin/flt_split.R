#!/usr/bin/env Rscript
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")

if (!requireNamespace("SIAMCAT", quietly = TRUE))
    BiocManager::install("SIAMCAT", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
    library(optparse); library(SIAMCAT)
})

opt_list <- list(
    make_option("--base", type="character"),
    make_option("--norm", type="character"),
    make_option("--cutoff", type="double"),
    make_option("--out", type="character"),
    make_option("--assoc_pdf", type="character"),
    make_option("--conf_pdf", type="character")
)
opt <- parse_args(OptionParser(option_list = opt_list))

sc <- readRDS(opt$base)

# Filter + association + confounders
sc <- filter.features(sc, filter.method='abundance', cutoff=opt$cutoff)
sc <- check.associations(sc, log.n0=1e-06, alpha=0.05)

# Plotting association features
pdf(opt$assoc_pdf, width=10, height=8)
association.plot(sc, max.show=50, sort.by='fc', panels=c('fc','prevalence','auroc'))
dev.off()

# Plotting confounders
# We let the function handle the PDF device opening/closing by passing the filename.
check.confounders(sc, meta.in=NULL, feature.type='filtered', fn.plot=opt$conf_pdf)


# Normalization params mapping (mirrors your R script)
norm_map <- list(
    'log.unit'  = list(norm.method='log.unit',  norm.param=list(log.n0=1e-05, n.p=1, norm.margin=1)),
    'log.std'   = list(norm.method='log.std',   norm.param=list(log.n0=1e-05, sd.min.q=.1)),
    'log.clr'   = list(norm.method='log.clr',   norm.param=list(log.n0=1e-06)),
    'rank.unit' = list(norm.method='rank.unit', norm.param=NULL),
    'rank.std'  = list(norm.method='rank.std',  norm.param=list(sd.min.q=0.1)),
    'std'       = list(norm.method='std',       norm.param=list(sd.min.q=0.1)),
    'pass'      = list(norm.method='pass',      norm.param=NULL)
)

nm <- norm_map[[opt$norm]]
if (is.null(nm)) stop("Unknown normalization: ", opt$norm)

if (is.null(nm$norm.param)) {
    sc <- normalize.features(sc, norm.method = nm$norm.method)
} else {
    sc <- normalize.features(sc, norm.method = nm$norm.method, norm.param = nm$norm.param)
}

sc <- create.data.split(sc, num.folds=5, num.resample=2)
saveRDS(sc, file=opt$out)


#!/usr/bin/env Rscript
# Score somatic SNV candidates with a trained SomaticXGB (XGBoost) model.
#
# Usage:
#   Rscript SomaticXGB.R <sample_readsName> <posfreq_file|-> <model_path>
#
# Run from the same analysis working directory used for pipeline 2 ("features").
# Expects:
#   mutect_a_vraw/<name>.raw.vcf.gz
#   varscan/<name>.raw.snp.sort.vcf.gz
#   varscan_novo/<name>.raw.snp.vcf.rc.tsv.gz
#
# posfreq_file: PosMix spike list for optional labels; pass "-" to skip labeling.
#
# Writes:
#   pred/<name>_pred.txt

library(data.table)
library(xgboost)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript SomaticXGB.R <name> <posfreq_file|-> <model_path>")
}
name <- args[1]
posfreq_file <- args[2]
model_path <- args[3]

if (identical(posfreq_file, "-") || !nzchar(posfreq_file)) {
  posfreq_file <- NULL
} else if (!file.exists(posfreq_file)) {
  stop("posfreq_file not found: ", posfreq_file)
}

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path)
}

ca <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", ca, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

source(file.path(script_dir, "SomaticXGB_features.R"))

dtm <- extract_somaticxgb_features(name, posfreq_file)

model <- xgb.load(model_path)
dtm$pred <- predict(model, data.matrix(dtm[, SOMATICXGB_FEATURES, with = FALSE]))

dir.create("pred", showWarnings = FALSE, recursive = TRUE)
fwrite(dtm, paste0("pred/", name, "_pred.txt"), sep = "\t")

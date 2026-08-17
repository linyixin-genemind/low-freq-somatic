#!/usr/bin/env Rscript
# Score somatic SNV candidates with a trained SomaticXGB (XGBoost) model.
#
# Usage:
#   Rscript SomaticXGB.R <sample_readsName> <posfreq|-> <model_path>
#
# Run from the pair's simulation working directory (pipeline-2 cwd for that pair),
# or let somaticxgb_predict.smk cd there via analysis_dir.
# Expects:
#   mutect_a_vraw/<name>.raw.vcf.gz
#   varscan/<name>.raw.snp.sort.vcf.gz
#   varscan_novo/<name>.raw.snp.vcf.rc.tsv.gz
#
# posfreq: PosMix spike list for optional labels; pass "-" to skip labeling.
#
# Writes:
#   pred/<name>_pred.txt

library(data.table)
library(xgboost)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript SomaticXGB.R <name> <posfreq|-> <model_path>")
}
name <- args[1]
posfreq <- args[2]
model_path <- args[3]

if (identical(posfreq, "-") || !nzchar(posfreq)) {
  posfreq <- NULL
} else if (!file.exists(posfreq)) {
  stop("posfreq not found: ", posfreq)
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

# Reuse the canonical feature extractor from step 3 (training).
features_r <- file.path(script_dir, "..", "3.SomaticXGB-training", "SomaticXGB_features.R")
if (!file.exists(features_r)) {
  stop("Cannot find shared extractor: ", normalizePath(features_r, mustWork = FALSE))
}
source(features_r)

dtm <- extract_somaticxgb_features(name, posfreq)

model <- xgb.load(model_path)
dtm$pred <- predict(model, data.matrix(dtm[, SOMATICXGB_FEATURES, with = FALSE]))

dir.create("pred", showWarnings = FALSE, recursive = TRUE)
fwrite(dtm, paste0("pred/", name, "_pred.txt"), sep = "\t")

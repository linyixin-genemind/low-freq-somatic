#!/usr/bin/env Rscript
# Score somatic SNV candidates with a trained SomaticXGB (XGBoost) model.
#
# Usage:
#   Rscript SomaticXGB.R <sample_readsName> <posfreq_file> [model_path]
#
# Expects, relative to the working directory:
#   mutect_a_vraw/<name>.raw.vcf.gz
#   varscan/<name>.raw.snp.sort.vcf.gz
#   varscan_novo/<name>.raw.snp.vcf.rc.tsv.gz
#
# Writes:
#   pred/<name>_pred.txt

library(data.table)
library(magrittr)
library(stringr)
library(tidyr)
library(tibble)
library(xgboost)

suppressPackageStartupMessages({
  # data.vcf provides read_vars(); install from the analysis environment if needed.
  library(data.vcf)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript SomaticXGB.R <name> <posfreq_file> [model_path]")
}
name <- args[1]
posfreq_file <- args[2]
model_path <- if (length(args) >= 3) args[3] else "SomaticXGB"

# Feature extraction for VarScan2-raw alleles force-called by Mutect2.
f <- function(name) {
  a <- read_vars(paste0("mutect_a_vraw/", name, ".raw.vcf.gz")) %>%
    .[grep("DP=", INFO)] %>%
    .[, TLOD := gsub(".*TLOD=", "", INFO)] %>%
    .[, INFO := gsub("[A-Z]", "", gsub("=", "", gsub(".*DP=", "", substr(INFO, 1, regexpr(";POPAF", INFO) - 1))))]

  dt <- a[, lapply(tstrsplit(INFO, ";", fixed = TRUE), as.character)] %>%
    .[, MBQ_REF := tstrsplit(V3, ",")[[1]]] %>%
    .[, MBQ_ALT := sub("^[^,]*,", "", V3)] %>%
    .[, MFRL_REF := tstrsplit(V4, ",")[[1]]] %>%
    .[, MFRL_ALT := sub("^[^,]*,", "", V4)] %>%
    .[, MMQ_REF := tstrsplit(V5, ",")[[1]]] %>%
    .[, MMQ_ALT := sub("^[^,]*,", "", V5)]

  dt <- cbind(dt[, !c("V3", "V4", "V5")], a[, .(TLOD, CHROM, POS, REF, ALT)])
  setnames(dt, old = names(dt)[1:5], new = c("DP", "ECNT", "MPOS", "NALOD", "NLOD"))

  a <- dt %>%
    as_tibble() %>%
    separate_rows(MPOS, NALOD, NLOD, MBQ_ALT, MFRL_ALT, MMQ_ALT, TLOD, ALT, sep = ",") %>%
    data.table() %>%
    .[, c_p := paste(CHROM, POS, sep = "_")]

  a <- a[nchar(REF) == nchar(ALT) & substr(REF, 1, 1) != substr(ALT, 1, 1) & substring(REF, 2) == substring(ALT, 2)] %>%
    .[, var := paste(CHROM, POS, substr(REF, 1, 1), substr(ALT, 1, 1), sep = ":")]

  v <- read_vars(paste0("varscan/", name, ".raw.snp.sort.vcf.gz")) %>%
    .[, var := paste(CHROM, POS, REF, ALT, sep = ":")] %>%
    .[, c_p := paste(CHROM, POS, sep = "_")]

  # Keep only alleles called by VarScan2.
  dt <- a[c_p %in% v$c_p] %>% .[var %in% v$var]

  pf <- fread(posfreq_file)
  if ("V4" %in% names(pf)) {
    pf[, var := paste(V1, V2, V4, V5, sep = ":")]
    dt[, label := as.numeric(var %in% pf$var)]
  }

  tsv <- fread(paste0("varscan_novo/", name, ".raw.snp.vcf.rc.tsv.gz")) %>%
    setnames(strsplit("CHROM:POS:REF:dp:b5:b3:base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end", ":")[[1]]) %>%
    .[, var := paste(CHROM, POS, REF, base, sep = ":")]

  sc <- c(
    "var", "avg_mapping_quality", "avg_basequality", "avg_pos_as_fraction",
    "avg_num_mismatches_as_fraction", "avg_sum_mismatch_qualities"
  )
  dtm <- merge(dt, tsv[, ..sc], by = "var")

  dtm[, id := gsub("_reads.*", "", name)]
  dtm[, reads := as.integer(gsub("M", "", gsub(".*_reads", "", name)))]
  dtm
}

dtm <- f(name)

features <- c(
  "DP", "ECNT", "MBQ_REF", "MBQ_ALT", "MFRL_REF", "MFRL_ALT", "MMQ_REF", "MMQ_ALT",
  "MPOS", "NALOD", "NLOD", "TLOD",
  "avg_mapping_quality", "avg_basequality", "avg_pos_as_fraction",
  "avg_num_mismatches_as_fraction", "avg_sum_mismatch_qualities"
)

for (col in features) {
  dtm[, (col) := as.numeric(get(col))]
}

if (!file.exists(model_path)) {
  stop("Model file not found: ", model_path)
}

model <- xgb.load(model_path)
dtm$pred <- predict(model, data.matrix(dtm[, features, with = FALSE]))

dir.create("pred", showWarnings = FALSE, recursive = TRUE)
fwrite(dtm, paste0("pred/", name, "_pred.txt"), sep = "\t")

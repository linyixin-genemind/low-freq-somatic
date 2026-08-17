# Intended use: with the working directory set to one pipeline-2 pair directory,
# call extract_somaticxgb_features(name, posfreq) to build and **return** one
# data.table for that job (does not write a file by itself).
#
# Callers:
#   - build_feature_matrix.R  — setwd(pair_dir) per pair, rbind returned tables
#     (used by train_SomaticXGB.ipynb and predict_SomaticXGB.ipynb)
# Do not duplicate this file under step 4.
#
# Inputs (paths relative to the current working directory = that pair dir):
#   mutect_a_vraw/<name>.raw.vcf.gz
#   varscan/<name>.raw.snp.sort.vcf.gz
#   varscan_novo/<name>.raw.snp.vcf.rc.tsv.gz
#   <posfreq> (optional)  — PosMix spiked sites → label column
#
# Returns: data.table of allele-level features (+ label when posfreq is given).

SOMATICXGB_FEATURES <- c(
  "DP", "ECNT", "MBQ_REF", "MBQ_ALT", "MFRL_REF", "MFRL_ALT", "MMQ_REF", "MMQ_ALT",
  "MPOS", "NALOD", "NLOD", "TLOD",
  "avg_mapping_quality", "avg_basequality", "avg_pos_as_fraction",
  "avg_num_mismatches_as_fraction", "avg_sum_mismatch_qualities"
)

.somaticxgb_require_inputs <- function(name) {
  needed <- c(
    paste0("mutect_a_vraw/", name, ".raw.vcf.gz"),
    paste0("varscan/", name, ".raw.snp.sort.vcf.gz"),
    paste0("varscan_novo/", name, ".raw.snp.vcf.rc.tsv.gz")
  )
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    stop(
      "Missing force-calling / feature inputs for '", name, "'.\n",
      "Expected under the current working directory:\n  ",
      paste(missing, collapse = "\n  "),
      "\nRun 2.PosMix-somatic-calling-pipeline with callers including \"features\" first."
    )
  }
}

extract_somaticxgb_features <- function(name, posfreq = NULL) {
  suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(stringr)
    library(tidyr)
    library(tibble)
    library(data.vcf)
  })

  .somaticxgb_require_inputs(name)

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

  if (!is.null(posfreq) && nzchar(posfreq) && posfreq != "-" && file.exists(posfreq)) {
    pf <- fread(posfreq)
    if ("V4" %in% names(pf)) {
      pf[, var := paste(V1, V2, V4, V5, sep = ":")]
    } else if (all(c("CHROM", "POS", "REF", "ALT") %in% names(pf))) {
      pf[, var := paste(CHROM, POS, REF, ALT, sep = ":")]
    } else {
      stop("Cannot parse variant keys from posfreq: ", posfreq)
    }
    dt[, label := as.numeric(var %in% pf$var)]
  } else {
    dt[, label := NA_real_]
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

  for (col in SOMATICXGB_FEATURES) {
    dtm[, (col) := as.numeric(get(col))]
  }
  dtm
}

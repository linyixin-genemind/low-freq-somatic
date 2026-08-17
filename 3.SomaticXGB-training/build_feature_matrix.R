#!/usr/bin/env Rscript
# Build the SomaticXGB training feature matrix from pipeline-2 force-calling outputs.
#
# Usage (run from any directory; paths in the config may be relative to this script
# or absolute — analysis_dir is where mutect_a_vraw/ varscan/ varscan_novo/ live):
#
#   Rscript build_feature_matrix.R config.training.json [output.tsv.gz]
#
# Each job must already have been produced by 2.PosMix-somatic-calling-pipeline
# with "features" in callers (HG001-003 training pairs).

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript build_feature_matrix.R <config.training.json> [output.tsv.gz]")
}

config_path <- normalizePath(args[1])
out_path <- if (length(args) >= 2) args[2] else "feature_matrix.tsv.gz"

cfg <- fromJSON(config_path, simplifyVector = TRUE)
analysis_dir <- cfg$analysis_dir
if (is.null(analysis_dir) || !dir.exists(analysis_dir)) {
  stop("config$analysis_dir must exist (directory with mutect_a_vraw/, varscan/, …): ", analysis_dir)
}

# Resolve feature extractor relative to this repo layout.
script_dir <- {
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
}
repo_root <- normalizePath(file.path(script_dir, ".."))
features_r <- file.path(repo_root, "4.SomaticXGB-predict", "SomaticXGB_features.R")
if (!file.exists(features_r)) {
  stop("Cannot find SomaticXGB_features.R at ", features_r)
}
source(features_r)

# Expand compact pairs × depths into jobs if "jobs" is absent.
expand_jobs <- function(cfg) {
  if (!is.null(cfg$jobs) && length(cfg$jobs)) {
    jobs <- as.data.table(cfg$jobs)
    # jsonlite may return a single job as a named list
    if (!is.data.frame(cfg$jobs) && is.null(jobs$name) && !is.null(cfg$jobs$name)) {
      jobs <- as.data.table(lapply(cfg$jobs, function(x) x))
    }
    return(jobs)
  }
  if (is.null(cfg$pairs) || is.null(cfg$tumor_reads)) {
    stop("Config must provide either 'jobs' or both 'pairs' and 'tumor_reads'")
  }
  pairs <- as.data.table(cfg$pairs)
  if (!nrow(pairs) && !is.null(cfg$pairs$group)) {
    pairs <- as.data.table(cfg$pairs)
  }
  reads <- cfg$tumor_reads
  pf_tmpl <- cfg$posfreq_template
  if (is.null(pf_tmpl)) {
    pf_tmpl <- "{pos_files}/train/{group}/pf_{id}_reads{reads}M.txt"
  }
  pos_files <- cfg$pos_files
  if (is.null(pos_files)) {
    pos_files <- file.path(repo_root, "pos_files")
  }

  jobs <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    p <- pairs[i]
    rbindlist(lapply(reads, function(r) {
      name <- paste0(p$id, "_reads", r, "M")
      pf <- pf_tmpl
      pf <- gsub("\\{pos_files\\}", pos_files, pf)
      pf <- gsub("\\{group\\}", p$group, pf)
      pf <- gsub("\\{id\\}", p$id, pf)
      pf <- gsub("\\{reads\\}", as.character(r), pf)
      data.table(
        name = name,
        posfreq = pf,
        group = as.character(p$group)
      )
    }))
  }))
  jobs
}

jobs <- expand_jobs(cfg)
if (!("posfreq" %in% names(jobs)) && ("posfreq_file" %in% names(jobs))) {
  message('WARNING: job column "posfreq_file" is deprecated; use "posfreq".')
  setnames(jobs, "posfreq_file", "posfreq")
}
required <- c("name", "posfreq")
miss <- setdiff(required, names(jobs))
if (length(miss)) stop("Each job needs columns: ", paste(required, collapse = ", "))

old_wd <- getwd()
setwd(analysis_dir)
on.exit(setwd(old_wd), add = TRUE)

parts <- list()
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i]
  name <- job$name
  pf <- job$posfreq
  if (!grepl("^/", pf)) {
    # Allow paths relative to repo root or config file directory.
    cand <- c(
      pf,
      file.path(dirname(config_path), pf),
      file.path(repo_root, pf)
    )
    hit <- cand[file.exists(cand)]
    if (!length(hit)) stop("posfreq not found for ", name, ": ", pf)
    pf <- normalizePath(hit[[1]])
  }
  message(sprintf("[%d/%d] %s", i, nrow(jobs), name))
  needed <- c(
    paste0("mutect_a_vraw/", name, ".raw.vcf.gz"),
    paste0("varscan/", name, ".raw.snp.sort.vcf.gz"),
    paste0("varscan_novo/", name, ".raw.snp.vcf.rc.tsv.gz")
  )
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    stop(
      "Missing force-calling / feature inputs under ", analysis_dir, ":\n  ",
      paste(missing, collapse = "\n  "),
      "\nRun 2.PosMix-somatic-calling-pipeline with callers including \"features\" first."
    )
  }
  dt <- extract_somaticxgb_features(name, pf)
  if ("group" %in% names(job) && !is.na(job$group)) dt[, group := as.character(job$group)]
  parts[[length(parts) + 1L]] <- dt
}

dtm <- rbindlist(parts, use.names = TRUE, fill = TRUE)
message("Feature matrix: ", nrow(dtm), " rows, ", sum(dtm$label == 1), " positives")

if (!grepl("^/", out_path)) {
  out_path <- file.path(old_wd, out_path)
}
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(dtm, out_path, sep = "\t")
message("Wrote ", out_path)

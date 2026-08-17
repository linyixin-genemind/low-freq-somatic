#!/usr/bin/env Rscript
# Build the SomaticXGB feature matrix from pipeline-2 force-calling outputs
# (shared by train_SomaticXGB.ipynb and predict_SomaticXGB.ipynb).
#
# No separate training config: each PosMix pair keeps the same config.json used for
# step 2 in that pair's working directory.
#
# Usage:
#   Rscript build_feature_matrix.R <pair_dir> [pair_dir ...] [output.tsv.gz]
#
# Example:
#   Rscript build_feature_matrix.R \
#     /path/to/train/12 \
#     /path/to/train/13 \
#     …/train/32 \
#     feature_matrix.tsv.gz
#
# Each <pair_dir> must contain:
#   config.json                          # the pipeline-2 config for that pair
#   mutect_a_vraw/{id}_reads{N}M.raw.vcf.gz
#   varscan/{id}_reads{N}M.raw.snp.sort.vcf.gz
#   varscan_novo/{id}_reads{N}M.raw.snp.vcf.rc.tsv.gz
#
# posfreq / tumor_reads / data[].id are read from that config.json (same keys as step 2).

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)


if (length(args) < 1) {
  stop(
    "Usage: Rscript build_feature_matrix.R <pair_dir> [pair_dir ...] [output.tsv.gz]\n",
    "Each pair_dir is a pipeline-2 working directory with its config.json."
  )
}

out_path <- "feature_matrix.tsv.gz"
pair_dirs <- args
if (grepl("\\.(tsv|txt|csv)(\\.gz)?$", args[length(args)], ignore.case = TRUE)) {
  out_path <- args[length(args)]
  pair_dirs <- args[-length(args)]
}
if (!length(pair_dirs)) {
  stop("Provide at least one pair directory (pipeline-2 cwd with config.json).")
}

# Expect to be invoked as:
#   Rscript /path/to/3.SomaticXGB-training/build_feature_matrix.R ...
# so sibling SomaticXGB_features.R is found even if cwd is step 4.
ca <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", ca, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
features_r <- file.path(script_dir, "SomaticXGB_features.R")
if (!file.exists(features_r)) {
  stop("SomaticXGB_features.R not found next to build_feature_matrix.R at ", script_dir)
}
source(features_r)

find_pair_config <- function(pair_dir) {
  config_path <- file.path(pair_dir, "config.json")
  if (!file.exists(config_path)) {
    stop(
      "Missing config.json under ", pair_dir,
      "\nEach pair directory must contain the pipeline-2 config.json."
    )
  }
  normalizePath(config_path)
}

# Absolute path as-is; relative path is resolved against the pair directory.
resolve_path <- function(path, pair_dir) {
  if (is.null(path) || !nzchar(path) || identical(path, "-")) return(NULL)
  if (grepl("^/", path)) {
    if (file.exists(path) || dir.exists(path)) return(normalizePath(path))
    return(NULL)
  }
  full <- file.path(pair_dir, path)
  if (file.exists(full) || dir.exists(full)) return(normalizePath(full))
  NULL
}

resolve_posfreq_for_job <- function(posfreq_path, sample_id, reads) {
  if (is.null(posfreq_path)) {
    stop("config must set \"posfreq\" (file or directory), same as pipeline 2")
  }
  if (dir.exists(posfreq_path)) {
    path <- file.path(posfreq_path, sprintf("pf_%s_reads%sM.txt", sample_id, reads))
    if (!file.exists(path)) {
      stop("No posfreq file for id=", sample_id, " reads=", reads, "M under ", posfreq_path,
           "\nExpected: ", path)
    }
    return(normalizePath(path))
  }
  if (file.exists(posfreq_path)) return(normalizePath(posfreq_path))
  stop("posfreq path not found (file or directory): ", posfreq_path)
}

# Infer pair tag from directory name (…/train/12 → "12") for the group column.
pair_group_from_dir <- function(pair_dir) {
  basename(normalizePath(pair_dir))
}

jobs_from_pair_dir <- function(pair_dir) {
  pair_dir <- normalizePath(pair_dir)
  if (!dir.exists(pair_dir)) stop("pair directory not found: ", pair_dir)
  config_path <- find_pair_config(pair_dir)
  cfg <- fromJSON(config_path, simplifyVector = TRUE)

  if (is.null(cfg$callers) || !("features" %in% cfg$callers)) {
    warning(
      "Pipeline-2 config for ", pair_dir,
      " does not list \"features\" in callers; force-calling outputs may be missing."
    )
  }

  posfreq_raw <- cfg$posfreq
  if (is.null(posfreq_raw)) posfreq_raw <- cfg$posfreq_file
  posfreq_path <- resolve_path(posfreq_raw, pair_dir)
  if (is.null(posfreq_path)) {
    stop("Cannot resolve posfreq from ", config_path, " relative to ", pair_dir, ": ", posfreq_raw)
  }

  if (is.null(cfg$tumor_reads) || !length(cfg$tumor_reads)) {
    stop("config tumor_reads missing in ", config_path)
  }
  if (is.null(cfg$data) || !length(cfg$data)) {
    stop("config data[] missing in ", config_path)
  }

  data <- cfg$data
  # jsonlite: one sample → named list; several → data.frame
  if (is.data.frame(data)) {
    ids <- as.character(data$id)
  } else if (!is.null(data$id)) {
    ids <- as.character(data$id)
  } else {
    ids <- vapply(data, function(x) as.character(x$id), character(1))
  }

  group <- pair_group_from_dir(pair_dir)
  reads <- cfg$tumor_reads

  rbindlist(lapply(ids, function(sid) {
    rbindlist(lapply(reads, function(r) {
      data.table(
        name = paste0(sid, "_reads", r, "M"),
        posfreq = resolve_posfreq_for_job(posfreq_path, sid, r),
        group = group,
        analysis_dir = pair_dir,
        config_path = config_path
      )
    }))
  }))
}

jobs <- rbindlist(lapply(pair_dirs, jobs_from_pair_dir))
message("Jobs: ", nrow(jobs), " from ", length(pair_dirs), " pair director",
        if (length(pair_dirs) == 1) "y" else "ies")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

parts <- list()
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i]
  name <- job$name
  adir <- job$analysis_dir
  pf <- job$posfreq

  message(sprintf("[%d/%d] %s  (dir=%s)", i, nrow(jobs), name, adir))
  setwd(adir)

  needed <- c(
    paste0("mutect_a_vraw/", name, ".raw.vcf.gz"),
    paste0("varscan/", name, ".raw.snp.sort.vcf.gz"),
    paste0("varscan_novo/", name, ".raw.snp.vcf.rc.tsv.gz")
  )
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    stop(
      "Missing force-calling / feature inputs under ", adir, ":\n  ",
      paste(missing, collapse = "\n  "),
      "\nRun 2.PosMix-somatic-calling-pipeline here with callers including \"features\"."
    )
  }
  dt <- extract_somaticxgb_features(name, pf)
  dt[, group := as.character(job$group)]
  parts[[length(parts) + 1L]] <- dt
}

setwd(old_wd)
dtm <- rbindlist(parts, use.names = TRUE, fill = TRUE)
message("Feature matrix: ", nrow(dtm), " rows, ", sum(dtm$label == 1), " positives")

if (!grepl("^/", out_path)) {
  out_path <- file.path(old_wd, out_path)
}
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(dtm, out_path, sep = "\t")
message("Wrote ", out_path)
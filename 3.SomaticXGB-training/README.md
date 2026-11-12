# SomaticXGB training

Build a labeled feature matrix from HG001–003 **force-calling** outputs, tune hyperparameters with 5-fold CV, and train the final XGBoost model.

## Prerequisite (required)

Training does **not** start from ordinary Mutect2 / Strelka2 / VarScan2 call sets alone. There is **no** separate training config — step 3 reuses each pair’s **pipeline-2** `config.json`.

Run `../2.PosMix-somatic-calling-pipeline/` **once per PosMix pair**, with `"features"` in `callers`, e.g.:

```json
"callers": ["varscan", "varscan_novo", "features"]
```

Pipeline-2 example configs use a **platform name** as `data[].id` (e.g. `"SURFSeq_5000"`), not a pair tag. Output files are therefore named `{id}_reads{N}M…`. If two pairs (e.g. `12` and `13`) share the same platform `id` in one working directory, they **overwrite each other**.

**Use a separate simulation working directory per pair**, and keep that pair’s step-2 config inside it. The tree below is an example of directories you create. `train_SomaticXGB.ipynb` reads the paths you write into `pair_dirs`; it does not create them, and running it unchanged does not produce the matrix or the model.

```
/path/to/your/train/12/
    config.json              ← same file used as snakemake --configfile
    mutect_a_vraw/…
    varscan/…
    varscan_novo/…
/path/to/your/train/13/
    config.json
    …
```

From each `config.json`, step 3 reads `data[].id`, `tumor_reads`, and `posfreq` (file or directory — same rules as pipeline 2). Without `"features"`, force-calling is skipped and training cannot run.

## Workflow

```
per-pair dirs (pipeline 2 config.json + "features" outputs)
        │
        ▼
train_SomaticXGB.ipynb          ← this directory
   • calls build_feature_matrix.R
        └─ sources SomaticXGB_features.R   ← library; not run alone
   • writes feature_matrix.tsv.gz, cv_results.tsv, med_auc.tsv,
     importance.tsv, SomaticXGB_retrained  (all in this directory)
   • 5-fold CV grid search (PR-AUC)
   • final model fit + xgb.save
```

`SomaticXGB_features.R` is a **library**, not a CLI. `build_feature_matrix.R` (used by both train and predict notebooks) `source()`s it and calls `extract_somaticxgb_features()`.

The training feature matrix is not in this repository. The notebook builds `feature_matrix.tsv.gz` from the per-pair force-calling outputs of pipeline 2. A file already written by that builder can be read with `fread` instead of rebuilding.

## R package `data.vcf`

`SomaticXGB_features.R` calls `library(data.vcf)` and uses `read_vars()` to read the force-calling and VarScan VCFs. The installed package is in `library/data.vcf` (built with R 4.2.3 on Linux x86-64). It is not on CRAN. Load it from this directory before running the feature script or either notebook:

```r
.libPaths(c("3.SomaticXGB-training/library", .libPaths()))
```

Use a path that reaches `library/` from the current working directory. The folder passed to `.libPaths()` is the parent of `data.vcf`, not `data.vcf` itself. The same R session also needs `data.table` (≥ 1.12.2), `Rcpp`, and `stringr`. R must be the 4.2 series; this copy is a compiled package, not source.

## Train in the notebook

1. Open `train_SomaticXGB.ipynb` in this directory (`3.SomaticXGB-training/`).
2. Replace the example `pair_dirs` with the pipeline-2 pair working directories you created (each with `config.json` + force-calling outputs). The notebook stops if those examples are left in place.
3. Run all cells:

| Step | What it does |
|------|----------------|
| Build matrix | Calls `build_feature_matrix.R` on `pair_dirs` (reads each pair’s step-2 config, `cd`s there, extracts features) |
| CV tune | Grid over `max_depth ∈ {3,5,7,9}`, `eta ∈ {0.01,0.05,0.1,0.2}`, `nrounds ∈ {100,200,400}`; 5-fold; median PR-AUC |
| Final train | Fits on **all** rows; default hyperparameters match the manuscript (`7` / `0.1` / `400`) |
| Outputs | In this directory: `feature_matrix.tsv.gz`, `cv_results.tsv`, `med_auc.tsv`, `importance.tsv`, `SomaticXGB_retrained` |

A `feature_matrix.tsv.gz` already written by this builder can be read with `fread` instead of running the build cell again.

Reference tables from the original analysis (`med_auc_0820.xls`, `importance_0820.txt`) remain in this directory for comparison.

## Released model

The manuscript release lives at `../models/SomaticXGB`, trained with:

- `max.depth = 7`
- `eta = 0.1`
- `nrounds = 400`
- `objective = "binary:logistic"`

Features match Supplementary Table 4 / `importance_0820.txt`. Retraining via the notebook writes `SomaticXGB_retrained` here.

## Prediction (no retraining)

Open `../4.SomaticXGB-predict/predict_SomaticXGB.ipynb`. It calls this directory’s `build_feature_matrix.R` on the test pair dirs, then loads a model and writes `predictions.tsv.gz`.

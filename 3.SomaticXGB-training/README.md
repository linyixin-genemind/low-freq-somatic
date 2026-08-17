# SomaticXGB training

Build a labeled feature matrix from HG001-003 **force-calling** outputs, tune hyperparameters with 5-fold CV, and train the final XGBoost model.

## Prerequisite (required)

Training does **not** start from ordinary Mutect2 / Strelka2 / VarScan2 call sets alone.

1. Run `../2.PosMix-somatic-calling-pipeline/` on the **training** PosMix pairs (`pos_files/train/{12,13,21,23,31,32}/`) with `"features"` in `callers`, e.g.:

```json
"callers": ["varscan", "varscan_novo", "features"]
```

2. That writes, for each `{id}_reads{reads}M`:

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz` — Mutect2 force-called VCF
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

Without `"features"`, force-calling is skipped and the steps below cannot run.

**Note on sample IDs:** output names are `{id}_reads{reads}M`, and directory-mode posfreq files must be `pf_{id}_reads{reads}M.txt`. Keep the same `id` in pipeline 2 and here. If several HG pairs share one working directory, their `id` values must still be unique **and** match the corresponding `pf_<id>_…` filenames (or run each pair in its own analysis directory).

## Workflow

```
force-calling VCFs (pipeline 2 + "features")
        │
        ▼
build_feature_matrix.R   ← shared extractor in 4.SomaticXGB-predict/SomaticXGB_features.R
        │
        ▼
feature_matrix.tsv.gz
        │
        ▼
notebooks/train_SomaticXGB.ipynb
   • 5-fold CV grid search (PR-AUC)
   • final model fit + xgb.save
```

### 1. Configure jobs

```bash
cp config.training.example.json config.training.json
# edit analysis_dir, pairs (group / id), tumor_reads
```

Compact form expands `pairs` × `tumor_reads` using `posfreq_template`. You can instead supply an explicit `jobs` array of `{name, posfreq, group}`.

### 2. Build the matrix (optional standalone)

```bash
Rscript build_feature_matrix.R config.training.json feature_matrix.tsv.gz
```

The training notebook runs this for you if `feature_matrix.tsv.gz` is built via the notebook cell.

### 3. Train in the notebook

Open `notebooks/train_SomaticXGB.ipynb` and run all cells:

| Step | What it does |
|------|----------------|
| Build matrix | Calls `build_feature_matrix.R` on HG001-003 force-calling outputs |
| CV tune | Grid over `max_depth ∈ {3,5,7,9}`, `eta ∈ {0.01,0.05,0.1,0.2}`, `nrounds ∈ {100,200,400}`; 5-fold; median PR-AUC |
| Final train | Fits on **all** rows; default hyperparameters match the manuscript (`7` / `0.1` / `400`) |
| Outputs | `cv_results.tsv`, `med_auc.tsv`, `importance.tsv`, model under `../models/` |

Reference tables from the original analysis (`med_auc_0820.xls`, `importance_0820.txt`) remain in this directory for comparison.

## Released model

`../models/SomaticXGB` was trained with:

- `max.depth = 7`
- `eta = 0.1`
- `nrounds = 400`
- `objective = "binary:logistic"`

Features match Supplementary Table 4 / `importance_0820.txt`.

## Prediction (no retraining)

Skip this directory and use `../4.SomaticXGB-predict/` with the released or newly saved model. Prediction inputs must also come from pipeline 2 with `"features"`.

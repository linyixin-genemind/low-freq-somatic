# SomaticXGB training

Manuscript setup: HG001–HG003 simulated tumor/normal pairs for training; HG004–HG005 for testing.

## Prerequisite (required)

Training does **not** start from ordinary Mutect2/Strelka2/VarScan2 call sets alone.

You must first run `../2.PosMix-somatic-calling-pipeline/` on the **training** PosMix pairs with `"features"` in `callers`, for example:

```json
"callers": ["varscan", "varscan_novo", "features"]
```

That step performs Mutect2 **force-calling** and writes the files the feature matrix is built from, including:

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz` — force-called VCF
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

If `"features"` was omitted in step 2, force-calling was skipped and these notebooks / feature assembly will not have the required inputs.

Use `pos_files/train/…` for training pairs (see `../pos_files/README.md`).

## Notebooks

| File | Role |
|------|------|
| `notebooks/tune_SomaticXGB.ipynb` | 5-fold CV grid search over `max_depth`, `eta`, `nrounds` (produces `med_auc_0820.xls`) |
| `notebooks/train_SomaticXGB.ipynb` | Final model fit with selected hyperparameters and `xgb.save` |

Selected hyperparameters (from `med_auc_0820.xls` / training notebook):

- `max.depth = 7`
- `eta = 0.1`
- `nrounds = 400`
- `objective = "binary:logistic"`

The released model is saved in the repo as `../models/SomaticXGB` (originally `model_0820`).

## Feature matrix

Training used a feature table assembled from PosMix + **force-calling** outputs (`dtm_0819.txt` in the original analysis directory). Features match Supplementary Table 4 / `importance_0820.txt`.

To rebuild features after running pipeline 2 with `"features"`:

1. Reuse the feature-extraction logic in `../4.SomaticXGB-predict/SomaticXGB.R` across all train pairs / depths.
2. Concatenate into one training table.
3. Run the tune → train notebooks above.

For scoring with a finished model (no retraining), skip this directory and use `../4.SomaticXGB-predict/` instead.

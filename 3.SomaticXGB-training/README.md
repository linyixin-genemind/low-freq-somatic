# SomaticXGB training

Manuscript setup: HG001–HG003 simulated tumor/normal pairs for training; HG004–HG005 for testing.

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

Training used a feature table assembled from PosMix + calling outputs (`dtm_0819.txt` in the original analysis directory). Features match Supplementary Table 4 / `importance_0820.txt`.

To rebuild features after running `2.PosMix-somatic-calling-pipeline` with `"features"` in `callers`, reuse the feature-extraction logic in `../4.SomaticXGB-predict/SomaticXGB.R` across all train pairs and concatenate.

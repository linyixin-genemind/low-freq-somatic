# SomaticXGB prediction

Score candidates with a trained / released model.

**There is no separate predict feature pipeline.** This step reuses training’s matrix builder:

`../3.SomaticXGB-training/build_feature_matrix.R` → feature table → load model → `pred` column.

The test feature matrix is not in this repository. The notebook builds it from the force-calling outputs of pipeline 2.

Before that script runs, point R at the bundled `data.vcf` package (R 4.2.x):

```r
.libPaths(c("3.SomaticXGB-training/library", .libPaths()))
```

See `../3.SomaticXGB-training/README.md`.

Primary entry point: `predict_SomaticXGB.ipynb` in this directory.

## Prerequisite (required)

Run `../2.PosMix-somatic-calling-pipeline/` **once per pair** to score, with `"features"` in `callers`, and keep that pair’s `config.json` in the working directory:

```
/path/to/your/test/45/
    config.json
    mutect_a_vraw/…
    varscan/…
    varscan_novo/…
```

`predict_SomaticXGB.ipynb` reads the paths you write into `pair_dirs`. It does not create those directories.

## Workflow

```
per-pair dirs (pipeline-2 config.json + force-calling outputs)
        │
        ▼
predict_SomaticXGB.ipynb          ← this directory
   • calls ../3.…/build_feature_matrix.R
        └─ sources SomaticXGB_features.R
   • writes feature_matrix.tsv.gz
   • xgb.load(model) → dtm$pred
   • writes predictions.tsv.gz
```

## Run

1. Open `predict_SomaticXGB.ipynb` (cwd = `4.SomaticXGB-predict/`).
2. Set `pair_dirs` to the test / held-out pair directories.
3. Set `model_path` (default `../models/SomaticXGB`, or a retrained model under step 3).
4. Run all cells.

Outputs in this directory: `feature_matrix.tsv.gz`, `predictions.tsv.gz`.

## Train vs predict

| Goal | Path |
|------|------|
| Train / CV | `../3.SomaticXGB-training/train_SomaticXGB.ipynb` |
| Predict | This notebook (same `build_feature_matrix.R`, then score) |

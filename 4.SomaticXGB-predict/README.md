# SomaticXGB prediction

Score candidates with a trained / released model (`../models/SomaticXGB` by default).

This directory does **not** run PosMix, callers, or force-calling. It only:

1. Reads pipeline-2 **force-calling** outputs  
2. Builds the same feature table used in training (`SomaticXGB_features.R`)  
3. Applies the XGBoost model → `pred/{id}_reads{reads}M_pred.txt.gz`

## Prerequisite (required)

Run `../2.PosMix-somatic-calling-pipeline/` in the **same working directory**, with `"features"` enabled:

```json
"callers": ["varscan", "varscan_novo", "features"]
```

Required files per `{id}` / `{reads}`:

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz`
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

Caller-only VCFs (Mutect2 / Strelka2 / VarScan **without** `"features"`) are **not** enough.

`data[].id` and `tumor_reads` in the predict config must match the names used in that pipeline-2 run.

## Config

See `config.example.json`. Relevant keys:

| Key | Role |
|-----|------|
| `data[].id` | Sample id → looks for `{id}_reads{reads}M` under `mutect_a_vraw/` etc. |
| `tumor_reads` | Depths to score (same list as pipeline 2) |
| `model` | Path to XGBoost model (default: `../models/SomaticXGB`) |
| `posfreq_file` | Optional. PosMix spike list for a `label` column (evaluation). Use `"-"` or omit for blind scoring |
| `data[].posfreq_file` | Optional per-sample override of the global `posfreq_file` |

BAM paths are **not** used here (features already come from force-calling outputs).

## Run

From the analysis directory that contains `mutect_a_vraw/`, `varscan/`, `varscan_novo/`:

```bash
snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
  --configfile config.json --cores 8
```

Or call the R script for one sample:

```bash
Rscript /path/to/4.SomaticXGB-predict/SomaticXGB.R \
  SURFSeq_5000_reads100M \
  ../pos_files/test/45/pf45_low.txt \
  ../models/SomaticXGB
```

Pass `-` instead of the posfreq path to skip labels.

## Train vs predict

| Goal | Path |
|------|------|
| Score with released / your model | This directory |
| Rebuild features + CV + retrain | `../3.SomaticXGB-training/` (same force-calling inputs) |

Feature extraction is shared: `SomaticXGB_features.R` ← also used by `3.…/build_feature_matrix.R`.

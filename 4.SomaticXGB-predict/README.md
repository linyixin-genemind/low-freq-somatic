# SomaticXGB prediction

Standalone Snakemake rules that apply a trained model (default: `../models/SomaticXGB`) to candidate variants.

This directory does **not** run PosMix, callers, or force-calling. It only scores files that already exist.

## Prerequisite (required)

Prediction inputs must come from `../2.PosMix-somatic-calling-pipeline/` run in the **same working directory**, with `"features"` enabled so Mutect2 **force-calling** and feature tables were produced.

Example step-2 `callers`:

```json
"callers": ["varscan", "varscan_novo", "features"]
```

Required files per sample / depth:

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz` — force-called VCF
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

Caller-only outputs (e.g. Mutect2 `filter.vcf.gz` or Strelka2 SNVs **without** `"features"`) are **not** sufficient for this pipeline.

If you only want conventional callers and no model scores, stop after pipeline 2 and omit `"features"` there; do not run this directory.

## Train vs predict

| Goal | Path |
|------|------|
| Use the released model (or your own `.model` file) | This directory (`4.…`) |
| Retrain / tune hyperparameters | `../3.SomaticXGB-training/` (still needs the same force-calling inputs from step 2) |

## Run

```bash
snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
  --configfile config.json --cores 8
```

Outputs: `pred/{id}_reads{reads}M_pred.txt.gz` (column `pred` = model score).

Optional config keys:

- `model` — path to XGBoost model (default: `../models/SomaticXGB`)
- `somaticxgb_rscript` — path to `SomaticXGB.R`

See `config.example.json`. Keep `data` / `tumor_reads` / `posfreq_file` consistent with the PosMix + `"features"` run that produced the inputs.

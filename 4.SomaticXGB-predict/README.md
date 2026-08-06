# SomaticXGB prediction

Standalone Snakemake rule extracted from the original PosMix/calling pipeline.
Applies the trained model in `../models/SomaticXGB` to candidate variants.

## Prerequisites

Run `../2.PosMix-somatic-calling-pipeline` in the **same working directory** with `"features"` enabled so these inputs exist:

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz`
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

## Run

```bash
snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
  --configfile config.json --cores 8
```

Outputs: `pred/{id}_reads{reads}M_pred.txt.gz` (column `pred` = model score).

Optional config keys:

- `model` — path to XGBoost model (default: `../models/SomaticXGB`)
- `somaticxgb_rscript` — path to `SomaticXGB.R`

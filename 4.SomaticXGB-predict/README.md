# SomaticXGB prediction

Score candidates with a trained / released model (`../models/SomaticXGB` by default).

This directory does **not** run PosMix, callers, or force-calling. It only:

1. Reads pipeline-2 **force-calling** outputs from a **per-pair** simulation directory  
2. Builds features via `../3.SomaticXGB-training/SomaticXGB_features.R`  
3. Applies the XGBoost model → scores under that pair’s `pred/` (and a copy under the snakemake cwd `pred/`)

## Prerequisite (required)

Pipeline-2 example configs use a **platform name** as `data[].id` (e.g. `SURFSeq_5000`). Run each PosMix pair in its **own** working directory (e.g. `{analysis_root}/test/45`), with `"features"` enabled:

```json
"callers": ["varscan", "varscan_novo", "features"]
```

Required files under that pair directory:

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz`
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

Caller-only VCFs (without `"features"`) are **not** enough.

## Config

See `config.example.json`. Relevant keys:

| Key | Role |
|-----|------|
| `analysis_dir` | Pair simulation directory (default: current working directory) |
| `data[].group` | Pair tag (`45`, `54`, …); combined with `id` so multiple pairs can share a platform name |
| `data[].id` | Platform id — must match pipeline 2 |
| `data[].analysis_dir` | Optional per-sample override of `analysis_dir` |
| `tumor_reads` | Depths to score (same list as pipeline 2) |
| `model` | Path to XGBoost model (default: `../models/SomaticXGB`) |
| `posfreq` | Optional. PosMix spike list for a `label` column. Use `"-"` or omit for blind scoring |
| `data[].posfreq` | Optional per-sample override of the global `posfreq` |

BAM paths are **not** used here.

## Run

Point `analysis_dir` at the pair directory (or `cd` there and omit it):

```bash
snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
  --configfile config.json --cores 8
```

Or call the R script **from inside** that pair directory:

```bash
cd /path/to/simulations/test/45
Rscript /path/to/4.SomaticXGB-predict/SomaticXGB.R \
  SURFSeq_5000_reads100M \
  /path/to/pos_files/test/45/pf45_low.txt \
  /path/to/models/SomaticXGB
```

Pass `-` instead of the posfreq path to skip labels.

## Train vs predict

| Goal | Path |
|------|------|
| Score with released / your model | This directory (per-pair `analysis_dir`) |
| Rebuild features + CV + retrain | `../3.SomaticXGB-training/` (reads each train pair’s pipeline-2 `config.json`) |

Feature extraction is **not duplicated here**: `SomaticXGB.R` sources `../3.SomaticXGB-training/SomaticXGB_features.R`.

# PosMix + somatic calling (+ optional force-calling)

Snakemake pipeline that combines three pieces:

1. **PosMix** — position-specific VAF mixing of tumor into normal background
2. **Mainstream somatic callers** — Mutect2, Strelka2, VarScan2 (± fpfilter)
3. **Force-calling / feature extraction** — only when `"features"` is listed in `callers`

SomaticXGB **training** and **prediction** are **not** in this directory. They live in `../3.SomaticXGB-training/` and `../4.SomaticXGB-predict/`. Both require outputs produced by enabling `"features"` here.

## Usage modes (choose via `callers`)

What this pipeline does is controlled almost entirely by the `callers` array in the config. Pick one mode below.

### Mode A — Callers only (no model work)

Use when you only want PosMix + conventional callers, and will **not** train or score with SomaticXGB.

```json
"callers": ["mutect2", "strelka2", "varscan", "varscan_novo"]
```

Omit any caller you do not need. Do **not** include `"features"`.  
Example config: `config.callers-only.example.json`.

### Mode B — Prepare inputs for SomaticXGB (force-calling required)

Use when you will either **train** a model (`3.SomaticXGB-training`) or **predict** with an existing model (`4.SomaticXGB-predict`).

`"features"` turns on:

- Mutect2 **force-calling** on VarScan2 alleles (`mutect_a_vraw/…`)
- Sorted VarScan2 SNP VCF
- bam-readcount TSV used as model features

Recommended config (matches the manuscript feature path):

```json
"callers": ["varscan", "varscan_novo", "features"]
```

You may also add `mutect2` / `strelka2` if you want platform comparison in the same run.  
Example config: `config.example.json` (features-enabled).

Without `"features"`, force-calling is skipped and later training/prediction steps will have nothing to read.

### Mode C — Callers + features in one run

Combine Mode A and B, e.g.:

```json
"callers": ["mutect2", "strelka2", "varscan", "varscan_novo", "features"]
```

Useful for full platform evaluation **and** SomaticXGB inputs from the same PosMix BAMs.

## `callers` reference

| Value | What it produces |
|-------|------------------|
| `mutect2` | Mutect2 filtered VCF |
| `strelka2` | Strelka2 somatic SNVs (via Manta indel candidates) |
| `varscan` | VarScan2 raw SNP/indel VCFs |
| `varscan_novo` | VarScan2 + fpfilter pass/fail VCFs |
| `features` | Force-calling VCF + sorted alleles + bam-readcount TSV for SomaticXGB |

Key `"features"` outputs (per `{id}` / `{reads}`):

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz` — force-called VCF (**required** for training and prediction)
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

## Downstream paths after this pipeline

```
2.PosMix… (with "features")
    ├─► 3.SomaticXGB-training/   # rebuild feature matrix → tune / train
    └─► 4.SomaticXGB-predict/    # score candidates with models/SomaticXGB
```

- **Train a new model**: run this pipeline on **training** PosMix pairs (`pos_files/train/…`) with `"features"`, then follow `../3.SomaticXGB-training/README.md`. The training notebooks need force-calling / feature files; caller-only VCFs are not enough.
- **Predict with the released (or your) model**: run this pipeline on the samples to score with `"features"`, then run `../4.SomaticXGB-predict/` in the **same working directory**. Prediction also requires the force-calling outputs above.

## Run

From an analysis working directory (BAM paths and `posfreq_file` set in the config):

```bash
# Optional: biobase Singularity wrapper if you use the original analysis image
sh /path/to/2.PosMix-somatic-calling-pipeline/run.sh \
  snakemake -s /path/to/2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk \
  --configfile config.json --cores 120
```

Or without the wrapper:

```bash
snakemake -s /path/to/2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk \
  --configfile config.json --cores 120
```

Point `posfreq_file` at files under `../pos_files/`. Panel BEDs are under `bed/panel/`.

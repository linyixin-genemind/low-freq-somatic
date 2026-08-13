# low-freq-somatic

Reproducibility repository for:

> *Assessing and Optimizing Low-Frequency Somatic Mutation Detection: A Multi-Platform High-Throughput Sequencing Perspective*

Pipeline order matches the manuscript analysis:

```
FASTQ → (1) fastp + BWA-MEM
           → (2) PosMix + somatic calling [+ force-calling if training or using model]
                → (3) SomaticXGB training
                     → (4) SomaticXGB prediction
```

## Repository layout

| Path | Description |
|------|-------------|
| `1.fastp-BWA-pipeline/` | Read QC (`fastp`) and alignment (`BWA-MEM` to GRCh37) |
| `2.PosMix-somatic-calling-pipeline/` | PosMix + Mutect2 / Strelka2 / VarScan2; optional `"features"` enables force-calling for SomaticXGB |
| `3.SomaticXGB-training/` | R notebooks for hyperparameter search and final model training |
| `4.SomaticXGB-predict/` | Apply a trained / released SomaticXGB model |
| `pos_files/` | Position–frequency tables for PosMix (`train/` = HG001–003, `test/` = HG004–005) |
| `models/SomaticXGB` | Released XGBoost model |

## Which path do you need?

| Goal | What to run | Config notes |
|------|-------------|--------------|
| Align only | `1.fastp-BWA-pipeline/` | See that directory’s `config.example.json` |
| PosMix + callers, **no** model | `2.…` | `callers` **without** `"features"` — see `config.callers-only.example.json` |
| Train SomaticXGB | `2.…` then `3.…` | In step 2, include `"features"` (force-calling). Training reads those feature / force-called files, not caller-only VCFs |
| Predict with existing model | `2.…` then `4.…` | Same: step 2 must include `"features"`; run prediction in the **same** analysis directory |

Details and example `callers` lists: `2.PosMix-somatic-calling-pipeline/README.md`.

## Manuscript mapping

- **Training simulations**: HG001 / HG002 / HG003 pairs → `pos_files/train/{12,13,21,23,31,32}/`
- **Test simulations**: HG004 / HG005 pairs → `pos_files/test/{45,54}/`
- **Model**: trained with `max.depth=7`, `eta=0.1`, `nrounds=400` (`3.SomaticXGB-training/`)

External tool used for mixing: [PosMix](https://github.com/liqg/posmix).

## Typical workflow

1. Align libraries with `1.fastp-BWA-pipeline/align.smk`.
2. Simulate and call with `2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk`:
   - Callers only → omit `"features"`.
   - Any SomaticXGB work → include `"features"` and set `posfreq_file` under `pos_files/`.
3. Either retrain under `3.SomaticXGB-training/`, or use `models/SomaticXGB`.
4. Score candidates with `4.SomaticXGB-predict/somaticxgb_predict.smk` (requires force-calling outputs from step 2).

See each subdirectory README for config examples and run commands.

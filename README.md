# low-freq-somatic

Reproducibility repository for:

> *Assessing and Optimizing Low-Frequency Somatic Mutation Detection: A Multi-Platform High-Throughput Sequencing Perspective*

Pipeline order matches the manuscript analysis:

```
FASTQ → (1) fastp + BWA-MEM
           → (2) PosMix + somatic calling (+ feature extraction)
               → (3) SomaticXGB training
                   → (4) SomaticXGB prediction
```

## Repository layout

| Path | Description |
|------|-------------|
| `1.fastp-BWA-pipeline/` | Read QC (`fastp`) and alignment (`BWA-MEM` to GRCh37) |
| `2.PosMix-somatic-calling-pipeline/` | PosMix simulation + Mutect2 / Strelka2 / VarScan2; optional `"features"` for model inputs. **No** SomaticXGB scoring |
| `3.SomaticXGB-training/` | R notebooks for hyperparameter search and final model training |
| `4.SomaticXGB-predict/` | Apply the trained SomaticXGB model |
| `pos_files/` | Position–frequency tables for PosMix (`train/` = HG001–003, `test/` = HG004–005) |
| `models/SomaticXGB` | Released XGBoost model (originally `model_0820`) |

## Manuscript mapping

- **Training simulations**: HG001 / HG002 / HG003 pairs → `pos_files/train/{12,13,21,23,31,32}/`
- **Test simulations**: HG004 / HG005 pairs → `pos_files/test/{45,54}/`
- **Model**: trained with `max.depth=7`, `eta=0.1`, `nrounds=400` (`3.SomaticXGB-training/`)

External tool used for mixing: [PosMix](https://github.com/liqg/posmix).

## Typical workflow

1. Align libraries with `1.fastp-BWA-pipeline/align.smk`.
2. Simulate and call with `2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk`, setting `posfreq_file` to the matching file under `pos_files/` and including `"features"` when preparing model inputs.
3. (Optional) Retrain with notebooks under `3.SomaticXGB-training/`, or use `models/SomaticXGB`.
4. Score candidates with `4.SomaticXGB-predict/somaticxgb_predict.smk`.

See each subdirectory README for config examples and run commands.

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

## Workflow paths

Shared start: align with `1.fastp-BWA-pipeline/align.smk`, then PosMix + calling via `2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk` (set `posfreq_file` under `pos_files/`). After that, choose a path:

| Goal | What to run next | Config notes |
|------|------------------|--------------|
| PosMix + callers only (**no** model) | Stop after `2.…` | `callers` **without** `"features"` — see `config.callers-only.example.json` |
| Train SomaticXGB | `2.…` → `3.SomaticXGB-training/` | In step 2 include `"features"` (force-calling). Training needs those feature / force-called files, not caller-only VCFs |
| Predict with existing / released model | `2.…` → `4.SomaticXGB-predict/` | Same: step 2 must include `"features"`; run prediction in the **same** analysis directory as step 2 |

For `callers` lists and run commands, see each subdirectory README (especially `2.PosMix-somatic-calling-pipeline/README.md`).

## Training and test sample details

- **Training**: HG001 / HG002 / HG003 PosMix pairs → `pos_files/train/{12,13,21,23,31,32}/`
- **Test**: HG004 / HG005 PosMix pairs → `pos_files/test/{45,54}/`
- **Released model hyperparameters**: `max.depth=7`, `eta=0.1`, `nrounds=400` (`3.SomaticXGB-training/`)

External tool used for mixing: [PosMix](https://github.com/liqg/posmix).

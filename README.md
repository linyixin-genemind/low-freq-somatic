# low-freq-somatic

Reproduction repository for the bioRxiv preprint:

> *Assessing and Optimizing Low-Frequency Somatic Mutation Detection: A Multi-Platform High-Throughput Sequencing Perspective*  
> https://www.biorxiv.org/content/10.64898/2026.05.28.728367v2

**Raw FASTQ data** are deposited in CNGBdb under project [CNP0009527](https://db.cngb.org/data_resources/project/CNP0009527/) (DOI: [10.26036/CNP0009527](https://doi.org/10.26036/CNP0009527)).

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
| `3.SomaticXGB-training/` | Feature extractor + build matrix (per-pair dirs) → CV tune → train SomaticXGB |
| `4.SomaticXGB-predict/` | Score with a trained / released model (reads each pair’s simulation dir; sources step-3 features) |
| `pos_files/` | Position–frequency tables for PosMix (`train/` = HG001-003, `test/` = HG004-005) |
| `models/SomaticXGB` | Released XGBoost model |

## Workflow paths

Shared start: align with `1.fastp-BWA-pipeline/align.smk`, then PosMix + calling via `2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk` in a **per-pair working directory** (pipeline-2 `data[].id` is a platform name). Set `"posfreq"` to a **file** (benchmark) or a **directory** under `pos_files/` (training). After that, choose a path:

| Goal | What to run next | Config notes |
|------|------------------|--------------|
| PosMix + callers only (**no** model) | Stop after `2.…` | `callers` **without** `"features"` — see `config.callers-only.example.json` |
| Train SomaticXGB | `2.…` **per pair dir** → `3.SomaticXGB-training/` | Step 2: `"features"` + per-pair `config.json`. Step 3 reads those configs (no separate training config) |
| Predict with existing / released model | `2.…` **per pair dir** → `4.SomaticXGB-predict/` | Step 2: `"features"`; step 4: `analysis_dir` = that pair’s simulation directory |

For `callers` lists and run commands, see each subdirectory README (especially `2.PosMix-somatic-calling-pipeline/README.md`).

## Training and test sample details

Pair IDs under `pos_files/` encode **donor → background** GIAB samples (first digit = tumor/donor, second = normal/background): e.g. `12` = HG001 reads mixed into HG002; `45` = HG004 into HG005.

- **Training**: HG001-003 PosMix pairs → `pos_files/train/{12,13,21,23,31,32}/`
- **Test**: HG004-005 PosMix pairs → `pos_files/test/{45,54}/`
- **Released model hyperparameters**: `max.depth=7`, `eta=0.1`, `nrounds=400` (`3.SomaticXGB-training/`)

External tool used for mixing: [PosMix](https://github.com/liqg/posmix).

# PosMix + somatic calling (+ optional force-calling)

Snakemake pipeline that combines three pieces:

1. **PosMix** — position-specific VAF mixing of tumor into normal background
2. **Mainstream somatic callers** — Mutect2, Strelka2, VarScan2 (± fpfilter)
3. **Force-calling / feature extraction** — only when `"features"` is listed in `callers`

SomaticXGB **training** and **prediction** are **not** in this directory. They live in `../3.SomaticXGB-training/` and `../4.SomaticXGB-predict/`. Both require outputs produced by enabling `"features"` here.

**Working directory per pair:** example configs use a **platform name** as `data[].id` (e.g. `"SURFSeq_5000"`), so outputs look like `mutect_a_vraw/SURFSeq_5000_reads100M.raw.vcf.gz`. Run **each PosMix pair in its own cwd** (e.g. `…/train/12`, `…/train/13`, `…/test/45`) so different pairs that share a platform id do not overwrite each other. Steps 3–4 read force-calling outputs from those per-pair directories.

## `posfreq`: file vs directory (auto-detected)

Config key: `"posfreq"`.


| What you set      | Mode                     | Behavior                                                                          |
| ----------------- | ------------------------ | --------------------------------------------------------------------------------- |
| **A single file** | Benchmark                | The **same** position–frequency table is used for every `id` × `tumor_reads`      |
| **A directory**   | Training-style diversity | For each job, load `pf_{id}_reads{N}M.txt` from that directory (see naming below) |


Examples shipped here:

- Benchmark / test-style: `"posfreq": "../pos_files/test/45/pf45_low.txt"` (`config.callers-only.example.json`)
- Training pair `12`: `"posfreq": "../pos_files/train/12"` (`config.example.json`) → matches `pf_SURFSeq_5000_reads100M.txt`, etc.

Directory matching: `pf_{id}_reads{N}M.txt` where `{id}` is `data[].id` and `{N}` is each value in `tumor_reads`.

Point `"posfreq"` at the **pair folder** (e.g. `train/12`), not at `train/` itself. How to build these tables: `../pos_files/README.md`.

## Usage modes (choose via `callers`)

What this pipeline does is controlled almost entirely by the `callers` array in the config. Pick one mode below.

### Mode A — Callers only (no model work)

Use when you only want PosMix + conventional callers, and will **not** train or score with SomaticXGB.

```json
"callers": ["mutect2", "strelka2", "varscan", "varscan_novo"]
```

Omit any caller you do not need. Do **not** include `"features"`.  
Example config: `config.callers-only.example.json` (single-file `posfreq` for a fair benchmark).

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

For **training diversity**, set `"posfreq"` to a pair directory under `../pos_files/train/…`.  
For **held-out test / benchmark scoring**, a single shared file (e.g. `test/45/pf45_low.txt`) is appropriate.

Example config: `config.example.json`. `"genome_fasta"` must be the reference FASTA; the pipeline does not supply a default path.

Without `"features"`, force-calling is skipped and later training/prediction steps will have nothing to read.

### Mode C — Callers + features in one run

Combine Mode A and B, e.g.:

```json
"callers": ["mutect2", "strelka2", "varscan", "varscan_novo", "features"]
```

Useful for full platform/caller evaluation **and** SomaticXGB inputs from the same PosMix BAMs.

## `callers` reference


| Value          | What it produces                                                      |
| -------------- | --------------------------------------------------------------------- |
| `mutect2`      | Mutect2 filtered VCF                                                  |
| `strelka2`     | Strelka2 somatic SNVs (via Manta indel candidates)                    |
| `varscan`      | VarScan2 SNP/indel VCFs after `bin/somaticfilter.pl`                  |
| `varscan_novo` | VarScan2 + fpfilter pass/fail VCFs                                    |
| `features`     | Force-calling VCF + sorted alleles + bam-readcount TSV for SomaticXGB |


Key `"features"` outputs (per `{id}` / `{reads}`):

- `mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz` — force-called VCF (**required** for training and prediction)
- `varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz`
- `varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz`

`features` turns bam-readcount output into that TSV with `bin/reshapeRC`. The executable and the two libraries it loads (`libseqrt.so`, `libomp.so`) are in `bin/`. The Snakemake rule sets `LD_LIBRARY_PATH` to `bin/`, so a normal `snakemake -s …/somatic_snv_posmix.smk` run does not need `reshapeRC` on `PATH`. The binary is Linux x86-64. To call it outside Snakemake:

```bash
export LD_LIBRARY_PATH=/path/to/2.PosMix-somatic-calling-pipeline/bin:${LD_LIBRARY_PATH}
bin/reshapeRC -i readcount.txt -g genome.fasta
```

`-i` is the bam-readcount text (the pipeline decompresses `*.rc.gz` and passes that stream). `-g` is the reference FASTA used for the read counts.

`bin/somaticfilter.pl` is plain Perl (no extra modules). The VarScan rule runs it on the somatic VCF before fpfilter and before the SomaticXGB candidate list is taken. Pass the VCF path, or `-` to read the VCF from standard input:

```bash
perl bin/somaticfilter.pl - < varscan.snp.vcf
```

## Run

From a **per-pair** analysis working directory. Keep that pair’s pipeline-2 config here (e.g. as `config.json`) — step 3 reuses it:

```bash
# example: /path/to/simulations/train/12   or   .../test/45
snakemake -s /path/to/2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk \
  --configfile config.json --cores 120
```

Panel BEDs are under `bed/`.
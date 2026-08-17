# pos_files

Position-and-frequency tables used by PosMix to spike in synthetic somatic SNVs.

Pipeline 2 auto-detects how you point at these tables via config `"posfreq"`:

- **File** → benchmark (same table for all depths / samples)
- **Directory** → per-`id` × per-depth files named `pf_<id>_reads<N>M.txt`

See `../2.PosMix-somatic-calling-pipeline/README.md`.

## Layout

```
pos_files/
  train/          # HG001-003 pairs (manuscript training set; diverse tables)
    12/ 13/ 21/ 23/ 31/ 32/
      pf_<id>_reads<N>M.txt          # PosMix tables (sites + target VAFs); <id> = pipeline data[].id
      only*.txt                     # donor-private site lists (no mixing ratios / VAFs)
  test/           # HG004-005 pairs (manuscript test / benchmark-style tables)
    45/pf45_low.txt
    54/pf54_low.txt
```

Pair IDs encode donor/background reference samples, e.g. `12` = HG001 reads mixed into HG002 background; `45` = HG004 into HG005.

This repository already includes **all** PosMix tables used for the manuscript training runs under `train/`. Prefer those files when reproducing training.

- **Use for PosMix / config `"posfreq"`**: `train/.../pf_*.txt` or `test/.../pf*.txt` (include target frequencies).
- **Site list only (no VAF)**: `train/<pair>/only*.txt` — donor-private germline sites for that pair, without mixing proportions.

## Format

Tab-delimited columns (no header), typically:

`CHROM  POS  P  REF  ALT  GT  f  b  GT2  group  AF`

`P` / `AF` are the target allele frequencies for PosMix.

## How these tables are generated (conceptual)

No generation scripts are shipped here—only the finished tables and this description.

### 1. GIAB donor-private germline sites

When simulating with GIAB reference samples, spiked “somatic” sites are taken from **germline variants private to the donor relative to the background**:

1. Download high-confidence small-variant VCFs from the [GIAB release FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/).
2. For a pair such as `12` (HG001 → HG002), take sites present in the **donor** (HG001) call set and absent from the **background** (HG002), restricted to your capture panel / analysis BED as needed.
3. Those donor-private heterozygous (or otherwise usable) SNVs become the candidate spike list (`only*.txt` in the train folders are related inputs from that step).

Useful entry points (use the `latest` or a pinned NIST version, GRCh37 to match this repo):

| Sample | Release directory |
|--------|-------------------|
| HG001 (NA12878) | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest |
| HG002 | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest |
| HG003 | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest |
| HG004 | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest |
| HG005 (Chinese trio) | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/latest |

Parent index: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/

### 2. Benchmark mode — identical content

For **fair caller / platform comparison**, every depth (and usually every replicate) should spike the **same** sites at the **same** target VAFs. Use one shared table (e.g. `test/45/pf45_low.txt`) and set pipeline config `"posfreq"` to that **file**.

### 3. Training mode — diversity across id × depth

For **model training**, increase diversity by drawing **different** site/VAF assignments for each `id` and each tumor read depth. Those tables live as:

`train/<pair>/pf_<id>_reads<N>M.txt`

Set pipeline `"posfreq"` to the **pair directory** (e.g. `../pos_files/train/12`). The workflow then picks `pf_<id>_reads<N>M.txt` automatically for each job.

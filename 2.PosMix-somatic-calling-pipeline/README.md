# PosMix + somatic calling

Snakemake pipeline for:

1. PosMix simulation (position-specific VAF mixing)
2. Mainstream somatic callers (Mutect2, Strelka2, VarScan2 ± fpfilter)
3. Feature extraction for SomaticXGB (`"features"` caller)

SomaticXGB **prediction** is intentionally **not** included here. After training (or using the released model), run `../4.SomaticXGB-predict/`.

## Callers

| Value in `callers` | Outputs |
|--------------------|---------|
| `mutect2` | Mutect2 filtered VCF |
| `strelka2` | Strelka2 somatic SNVs (via Manta indel candidates) |
| `varscan` | VarScan2 raw SNP/indel VCFs |
| `varscan_novo` | VarScan2 + fpfilter pass/fail |
| `features` | Force-calling + bam-readcount TSV used by SomaticXGB |

For manuscript training-style runs use e.g. `["varscan", "varscan_novo", "features"]`.  
For full platform evaluation include Mutect2/Strelka2 as needed.

## Run

From an analysis working directory (BAM paths and `posfreq_file` in the config):

```bash
# With the biobase Singularity image (same as original analysis):
sh /path/to/2.PosMix-somatic-calling-pipeline/run.sh \
  snakemake -s /path/to/2.PosMix-somatic-calling-pipeline/somatic_snv_posmix.smk \
  --configfile config.json --cores 120
```

See `config.example.json`. Point `posfreq_file` at files under `../pos_files/`.
Panel BEDs are under `bed/panel/`.

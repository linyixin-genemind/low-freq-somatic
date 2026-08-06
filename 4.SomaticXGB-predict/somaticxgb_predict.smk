# Apply the trained SomaticXGB model to PosMix + calling outputs.
#
# Prerequisites: run 2.PosMix-somatic-calling-pipeline with callers that include
# "features" (and typically "varscan" / "varscan_novo") so that the following
# exist for each sample/depth:
#   mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz
#   varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz
#   varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz
#
# Usage (from the same analysis directory used for PosMix/calling):
#   snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
#       --configfile config.json -j <cores>

import os

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(workflow.snakefile), ".."))
PIPELINE_DIR = os.path.dirname(workflow.snakefile)

model_path = config.get("model", os.path.join(REPO_ROOT, "models", "SomaticXGB"))
rscript = config.get("somaticxgb_rscript", os.path.join(PIPELINE_DIR, "SomaticXGB.R"))
posfreq_file = config["posfreq_file"]

ls = []
for pair in config["data"]:
    for reads in config["tumor_reads"]:
        ls.append(f'pred/{pair["id"]}_reads{reads}M_pred.txt.gz')


rule all:
    input:
        ls


rule SomaticXGB:
    """Score VarScan2 candidates with the trained SomaticXGB model."""
    input:
        "mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz",
        "varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz",
        "varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz",
    output:
        "pred/{id}_reads{reads}M_pred.txt.gz",
    benchmark:
        "logs/{id}_reads{reads}M_somaticxgb.benchmark.txt"
    params:
        pf=posfreq_file,
        name="{id}_reads{reads}M",
        model=model_path,
        rscript=rscript,
    shell:
        """
        mkdir -p pred
        Rscript {params.rscript} {params.name} {params.pf} {params.model}
        bgzip -f pred/{params.name}_pred.txt
        """

# Apply a trained SomaticXGB model to PosMix + force-calling outputs.
#
# Required upstream: 2.PosMix-somatic-calling-pipeline with "features" in callers
# (run prediction from that same working directory). Without force-calling:
#   mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz
#   varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz
#   varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz
# will be missing and this pipeline cannot run.
#
# Usage:
#   snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
#       --configfile config.json -j <cores>

import os

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(workflow.snakefile), ".."))
PIPELINE_DIR = os.path.dirname(workflow.snakefile)


def _resolve(path, default=None):
    if path is None or path == "" or path == "-":
        return default
    if os.path.isabs(path):
        return path
    # Prefer path relative to cwd (analysis dir); fall back to repo root.
    cand = [
        os.path.abspath(path),
        os.path.abspath(os.path.join(REPO_ROOT, path)),
        os.path.abspath(os.path.join(PIPELINE_DIR, path)),
    ]
    for p in cand:
        if os.path.exists(p):
            return p
    return cand[0]


model_path = _resolve(
    config.get("model"),
    default=os.path.join(REPO_ROOT, "models", "SomaticXGB"),
)
rscript = _resolve(
    config.get("somaticxgb_rscript"),
    default=os.path.join(PIPELINE_DIR, "SomaticXGB.R"),
)
# Global default; each sample may override with data[].posfreq
# (legacy alias "posfreq_file" still accepted).
default_posfreq = config.get("posfreq")
if default_posfreq is None:
    default_posfreq = config.get("posfreq_file", "-")
    if default_posfreq != "-":
        print('[posfreq] WARNING: config key "posfreq_file" is deprecated; use "posfreq".')

sample_posfreq = {}
for pair in config["data"]:
    pf = pair.get("posfreq")
    if pf is None:
        pf = pair.get("posfreq_file", default_posfreq)
    if pf not in (None, "", "-"):
        pf = _resolve(pf)
    else:
        pf = "-"
    sample_posfreq[pair["id"]] = pf

ls = []
for pair in config["data"]:
    for reads in config["tumor_reads"]:
        ls.append(f'pred/{pair["id"]}_reads{reads}M_pred.txt.gz')


rule all:
    input:
        ls


rule SomaticXGB:
    """Extract features from force-calling outputs and score with SomaticXGB."""
    input:
        mutect="mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz",
        varscan="varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz",
        rc="varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz",
    output:
        "pred/{id}_reads{reads}M_pred.txt.gz",
    benchmark:
        "logs/{id}_reads{reads}M_somaticxgb.benchmark.txt"
    params:
        pf=lambda wildcards: sample_posfreq[wildcards.id],
        name="{id}_reads{reads}M",
        model=model_path,
        rscript=rscript,
    shell:
        """
        mkdir -p pred logs
        test -f {params.model} || {{ echo "Model not found: {params.model}" >&2; exit 1; }}
        Rscript {params.rscript} {params.name} {params.pf} {params.model}
        bgzip -f pred/{params.name}_pred.txt
        """

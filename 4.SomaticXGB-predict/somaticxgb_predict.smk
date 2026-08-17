# Apply a trained SomaticXGB model to PosMix + force-calling outputs.
#
# Required upstream: 2.PosMix-somatic-calling-pipeline with "features" in callers,
# in a **per-pair** working directory (pipeline-2 data[].id is the platform name).
#
# Typical usage — cd into that pair's simulation directory:
#   snakemake -s /path/to/4.SomaticXGB-predict/somaticxgb_predict.smk \
#       --configfile config.json -j <cores>
#
# Or set analysis_dir / data[].analysis_dir / data[].group so several pair
# directories can be scored in one run without colliding on platform ids.

import os

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(workflow.snakefile), ".."))
PIPELINE_DIR = os.path.dirname(workflow.snakefile)


def _resolve(path, default=None):
    if path is None or path == "" or path == "-":
        return default
    if os.path.isabs(path):
        return path
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

default_analysis_dir = config.get("analysis_dir")
if default_analysis_dir in (None, "", "-"):
    default_analysis_dir = os.getcwd()
else:
    default_analysis_dir = _resolve(default_analysis_dir)

default_posfreq = config.get("posfreq")
if default_posfreq is None:
    default_posfreq = config.get("posfreq_file", "-")
    if default_posfreq != "-":
        print('[posfreq] WARNING: config key "posfreq_file" is deprecated; use "posfreq".')

# key = "{group}__{id}" so the same platform id can appear for multiple pairs.
samples = {}
for pair in config["data"]:
    sid = pair["id"]
    group = str(pair.get("group", "run"))
    key = f"{group}__{sid}"
    adir = pair.get("analysis_dir", default_analysis_dir)
    if adir in (None, "", "-"):
        adir = default_analysis_dir
    else:
        adir = _resolve(adir)
    pf = pair.get("posfreq")
    if pf is None:
        pf = pair.get("posfreq_file", default_posfreq)
    if pf not in (None, "", "-"):
        pf = _resolve(pf)
    else:
        pf = "-"
    if key in samples:
        raise ValueError(f"Duplicate sample key {key}; set distinct data[].group per pair.")
    samples[key] = {
        "id": sid,
        "group": group,
        "analysis_dir": adir,
        "posfreq": pf,
    }

ls = [
    f"pred/{key}_reads{reads}M_pred.txt.gz"
    for key in samples
    for reads in config["tumor_reads"]
]


rule all:
    input:
        ls


rule SomaticXGB:
    """Score one platform×depth job; reads force-calling files from that pair's dir."""
    output:
        "pred/{key}_reads{reads}M_pred.txt.gz",
    wildcard_constraints:
        reads="[0-9]+",
        key="[^/]+",
    benchmark:
        "logs/{key}_reads{reads}M_somaticxgb.benchmark.txt"
    params:
        rscript=rscript,
        model=model_path,
        adir=lambda wildcards: samples[wildcards.key]["analysis_dir"],
        sid=lambda wildcards: samples[wildcards.key]["id"],
        pf=lambda wildcards: samples[wildcards.key]["posfreq"],
        name=lambda wildcards: f'{samples[wildcards.key]["id"]}_reads{wildcards.reads}M',
    input:
        mutect=lambda wildcards: os.path.join(
            samples[wildcards.key]["analysis_dir"],
            "mutect_a_vraw",
            f'{samples[wildcards.key]["id"]}_reads{wildcards.reads}M.raw.vcf.gz',
        ),
        varscan=lambda wildcards: os.path.join(
            samples[wildcards.key]["analysis_dir"],
            "varscan",
            f'{samples[wildcards.key]["id"]}_reads{wildcards.reads}M.raw.snp.sort.vcf.gz',
        ),
        rc=lambda wildcards: os.path.join(
            samples[wildcards.key]["analysis_dir"],
            "varscan_novo",
            f'{samples[wildcards.key]["id"]}_reads{wildcards.reads}M.raw.snp.vcf.rc.tsv.gz',
        ),
    shell:
        """
        mkdir -p pred logs "{params.adir}/pred"
        test -f {params.model} || {{ echo "Model not found: {params.model}" >&2; exit 1; }}
        (
          cd "{params.adir}"
          Rscript {params.rscript} {params.name} {params.pf} {params.model}
          bgzip -f "pred/{params.name}_pred.txt"
        )
        cp -f "{params.adir}/pred/{params.name}_pred.txt.gz" "{output}"
        """

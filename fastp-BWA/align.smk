# Usage: snakemake -s align.smk --configfile config.example.json -j <threads>

import os

genome_fasta = config["genome_fasta"]
bed_file = config.get("bed_file", "")
fastp_args = config.get("fastp", "")

fqs = {p["sampleid"]: p for p in config["data"]}
for sid, p in fqs.items():
    p["ispe"] = len(p.get("read2", [])) > 0
sampleids = list(fqs.keys())


rule all:
    input:
        expand("pipeline/bwamem/{sample}/{sample}.bam", sample=sampleids),
        expand("pipeline/bwamem/{sample}/{sample}.bam.bai", sample=sampleids),


rule catfq:
    """Merge multi-lane FASTQs (or symlink a single pair)."""
    input:
        fq1=lambda wildcards: fqs[wildcards.sample]["read1"],
        fq2=lambda wildcards: fqs[wildcards.sample].get("read2", []),
    output:
        fq1=temp("pipeline/reads/{sample}/{sample}_R1.fq.gz"),
        fq2=temp("pipeline/reads/{sample}/{sample}_R2.fq.gz"),
    run:
        r1, r2 = list(input.fq1), list(input.fq2)
        ispe = len(r2) > 0
        if len(r1) < 1:
            raise Exception(f"{wildcards.sample}: read1 is empty")
        if ispe and len(r1) != len(r2):
            raise Exception(f"{wildcards.sample}: read1/read2 count mismatch")
        if not ispe:
            shell("touch {output.fq2}")
        if len(r1) == 1:
            shell("ln -srf {r1[0]} {output.fq1}")
            if ispe:
                shell("ln -srf {r2[0]} {output.fq2}")
        else:
            shell("cat {input.fq1} > {output.fq1}")
            if ispe:
                shell("cat {input.fq2} > {output.fq2}")


rule fastp:
    input:
        fq1="pipeline/reads/{sample}/{sample}_R1.fq.gz",
        fq2="pipeline/reads/{sample}/{sample}_R2.fq.gz",
    output:
        fq1=temp("pipeline/fastp/{sample}/{sample}_R1.fq.gz"),
        fq2=temp("pipeline/fastp/{sample}/{sample}_R2.fq.gz"),
        json="pipeline/fastp/{sample}/{sample}.json",
        html="pipeline/fastp/{sample}/{sample}.html",
    threads: 6
    params:
        extra=lambda wildcards: (
            fastp_args + " --detect_adapter_for_pe"
            if fqs[wildcards.sample]["ispe"]
            else fastp_args
        ),
    shell:
        """
        if [ -s "{input.fq2}" ]; then
            fastp --dont_eval_duplication --thread {threads} \
                -i {input.fq1} -o {output.fq1} \
                -I {input.fq2} -O {output.fq2} \
                --json {output.json} --html {output.html} \
                {params.extra}
        else
            fastp --dont_eval_duplication --thread {threads} \
                -i {input.fq1} -o {output.fq1} \
                --json {output.json} --html {output.html} \
                {params.extra}
            touch {output.fq2}
        fi
        """


rule bwamem:
    input:
        fq1="pipeline/fastp/{sample}/{sample}_R1.fq.gz",
        fq2="pipeline/fastp/{sample}/{sample}_R2.fq.gz",
    output:
        bam="pipeline/bwamem/{sample}/{sample}.bam",
        bai="pipeline/bwamem/{sample}/{sample}.bam.bai",
    threads: 16
    params:
        ref=genome_fasta,
        rg=r'-R "@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"',
        # Matches manuscript Methods: BWA-MEM -K 100000000 -k 32 -M
        bwa_opts="-K 100000000 -k 32 -M",
    shell:
        """
        if [ -s "{input.fq2}" ]; then
            bwa mem -t {threads} {params.bwa_opts} {params.rg} {params.ref} \
                {input.fq1} {input.fq2} \
                | samtools sort -@ 4 -m 2G -O BAM -o {output.bam}
        else
            bwa mem -t {threads} {params.bwa_opts} {params.rg} {params.ref} \
                {input.fq1} \
                | samtools sort -@ 4 -m 2G -O BAM -o {output.bam}
        fi
        samtools index {output.bam}
        """

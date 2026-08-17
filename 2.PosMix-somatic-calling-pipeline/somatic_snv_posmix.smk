# PosMix + somatic callers + optional force-calling for SomaticXGB.
#
# Controlled by config["callers"]:
#   mutect2 / strelka2 / varscan / varscan_novo  — mainstream callers only
#   features — Mutect2 force-calling + bam-readcount TSV for training/prediction
#
# posfreq:
#   file      → benchmark (same table for every id × depth)
#   directory → pf_{id}_reads{N}M.txt per job (training diversity)
#
# Callers-only (no model): omit "features".
# Train (3.SomaticXGB-training) or predict (4.SomaticXGB-predict): include "features".
# See README.md for usage modes and example configs.

import os, random, json

os.system("mkdir -p posmix")

for pair in config["data"]:
    bam_normal = pair["normal"]
    bam_tumor = pair["tumor"]
    bam_id = pair["id"]
    
    N_raw_stats = f"{bam_normal}.flagstat.txt"
    if not os.path.exists(N_raw_stats):
        N_raw_stats = f"posmix/{bam_id}_N.raw.flagstat.txt"
        if not os.path.exists(N_raw_stats):
            os.system(f"samtools flagstat -@ 8 {bam_normal} > {N_raw_stats}")

    # Take the read count from the first line of the .flagstat.txt file
    # All read counts in these files are counted as single-end (SE)
    normal_total = int(os.popen(f"sed -n '1p' {N_raw_stats} | grep -o '^[0-9]*'").read())
    if normal_total < max(max(list(map(float, config["tumor_reads"]))), float(config["normal_reads"])) * 1000000:
        raise Exception(f"Insufficient reads in {bam_normal}, please check again!")

ls = []
seed = int(config["seed"])
N_pileup = {}
N_tidy = {}
bed_nsplit = 20
splitids = [str(i).zfill(4) for i in range(bed_nsplit)]
seed_dict = {}
bams = {}

for pair in config["data"]:
    bams[pair["id"]] = pair
    N_tidy[pair["id"]] = f'posmix/{pair["id"]}_N.tidy.bam'
    N_pileup[pair["id"]] = f'varscan/{pair["id"]}_N.pileup.gz'
    seed_dict[f'{pair["id"]}_N'] = seed
    seed += 2

for pair in config["data"]:
    for reads in config["tumor_reads"]:

        if "strelka2" in config["callers"]:
            # ls.append(f'manta/{pair["id"]}_reads{reads}M/results/variants/candidateSmallIndels.vcf.gz')
            ls.append(f'strelka/{pair["id"]}_reads{reads}M/results/variants/somatic.snvs.vcf.gz')

        if "mutect2" in config["callers"]:
            for s in splitids:
                ls.append(f'mutect/{pair["id"]}_reads{reads}M.tmpvcf/{s}.vcf.gz')
            ls.append(f'mutect/{pair["id"]}_reads{reads}M.raw.vcf.gz')
            ls.append(f'mutect/{pair["id"]}_reads{reads}M.raw.vcf.gz.stats')
            ls.append(f'mutect/{pair["id"]}_reads{reads}M.filter.vcf.gz')

        if "varscan" in config["callers"]:
            ls.append(f'varscan/{pair["id"]}_reads{reads}M.raw.snp.vcf.gz')
            ls.append(f'varscan/{pair["id"]}_reads{reads}M.raw.indel.vcf.gz')
        
        if "varscan_novo" in config["callers"]:
            ls.append(f'varscan_novo/{pair["id"]}_reads{reads}M.snp.pass.vcf.gz')
            ls.append(f'varscan_novo/{pair["id"]}_reads{reads}M.indel.pass.vcf.gz')
            ls.append(f'varscan_novo/{pair["id"]}_reads{reads}M.snp.failed.vcf.gz')
            ls.append(f'varscan_novo/{pair["id"]}_reads{reads}M.indel.failed.vcf.gz')

        # Enabled only when "features" is in callers. Required before
        # 3.SomaticXGB-training or 4.SomaticXGB-predict (force-calling + TSV).
        if "features" in config["callers"]:
            ls.append(f'mutect_a_vraw/{pair["id"]}_reads{reads}M.raw.vcf.gz')
            ls.append(f'varscan/{pair["id"]}_reads{reads}M.raw.snp.sort.vcf.gz')
            ls.append(f'varscan_novo/{pair["id"]}_reads{reads}M.raw.snp.vcf.rc.tsv.gz')

        seed_dict[f'{pair["id"]}_T_reads{reads}M'] = seed
        seed += 2

ls.append('mutect/intervals.list')

jsObj = json.dumps(seed_dict)
fileObject = open('seed.json', 'w')
fileObject.write(jsObj)
fileObject.close()

# Prefer explicit genome_fasta; fall back to biobase layout used in the original analysis.
genome_dir = config.get(
    "genome_fasta",
    f"/biobase/data/genome/{config['genome']}/genome.fasta",
)

# ---------------------------------------------------------------------------
# Posfreq path: file = benchmark (same table for every id×depth);
#               directory = training-style (match pf_{id}_reads{N}M.txt).
# Config key: "posfreq" (legacy alias "posfreq_file" still accepted).
# ---------------------------------------------------------------------------
posfreq_path = config.get("posfreq")
if posfreq_path is None and "posfreq_file" in config:
    print('[posfreq] WARNING: config key "posfreq_file" is deprecated; use "posfreq".')
    posfreq_path = config["posfreq_file"]
if posfreq_path is None:
    raise Exception('Config must set "posfreq" to a file (benchmark) or a directory (per id×depth tables).')
posfreq_path = os.path.abspath(posfreq_path)

if os.path.isdir(posfreq_path):
    posfreq_mode = "directory"
elif os.path.isfile(posfreq_path):
    posfreq_mode = "file"
else:
    raise Exception(f"posfreq path not found (file or directory): {posfreq_path}")


def resolve_posfreq(sample_id, reads):
    """Return the PosMix position–frequency table for one id × depth."""
    if posfreq_mode == "file":
        return posfreq_path
    path = os.path.join(posfreq_path, f"pf_{sample_id}_reads{reads}M.txt")
    if os.path.isfile(path):
        return path
    raise Exception(
        f"No posfreq file for id={sample_id} reads={reads}M under {posfreq_path}.\n"
        f"Expected: {path}\n"
        "Directory mode requires names pf_<id>_reads<N>M.txt (see ../pos_files/README.md)."
    )


# Resolve and validate all tables needed for this config; collect for cp.txt.
posfreq_by_job = {}
all_posfreq_files = []
for pair in config["data"]:
    for reads in config["tumor_reads"]:
        pf = resolve_posfreq(pair["id"], reads)
        posfreq_by_job[(pair["id"], str(reads))] = pf
        all_posfreq_files.append(pf)
all_posfreq_files = sorted(set(all_posfreq_files))

print(f"[posfreq] mode={posfreq_mode} path={posfreq_path} ({len(all_posfreq_files)} table(s))")


rule all:
    input:
        ls

rule cp:
    input: all_posfreq_files
    output: "cp.txt"
    shell: "cut -f1,2 {input} | sort -u > {output}"

rule N_sample:
    input:
        N = lambda wildcards: bams[wildcards.id]['normal']
    output:
        temp("posmix/{id}_N.bam")
    threads: 4
    benchmark: "logs/{id}_N_sample.benchmark.txt"
    params:
        n = float(config["normal_reads"]) * 1000000
    run:
        N_raw_stats = f"{input.N}.flagstat.txt"
        if not os.path.exists(N_raw_stats):
            N_raw_stats = f"posmix/{wildcards.id}_N.raw.flagstat.txt"
        ntotal = int(os.popen(f"sed -n '1p' {N_raw_stats} | grep -o '^[0-9]*'").read())
        p = params.n / ntotal
        seed = seed_dict[f"{wildcards.id}_N"]
        shell("samtools view --no-PG -@ 8 -s " + str(seed + p) + " -b -o {output} {input.N}")
        shell("samtools index {output}")

rule posmix:
    input:
        posfreq = lambda wildcards: posfreq_by_job[(wildcards.id, str(wildcards.reads))],
        T = lambda wildcards: bams[wildcards.id]['tumor'],
        N = lambda wildcards: bams[wildcards.id]['normal'],
        cp = "cp.txt"
    output:
        receptor = temp("posmix/{id}_T_reads{reads}M.receptor.bam"),
        posmix_raw = temp("posmix/{id}_T_reads{reads}M.raw.bam"),
        posmix = "posmix/{id}_T_reads{reads}M.bam"
    threads: 8
    benchmark: "logs/{id}_reads{reads}M_posmix.benchmark.txt"
    run:
        N_raw_stats = f"{input.N}.flagstat.txt"
        if not os.path.exists(N_raw_stats):
            N_raw_stats = f"posmix/{wildcards.id}_N.raw.flagstat.txt"
        
        n = float(wildcards.reads) * 1000000            
        n_total_N = int(os.popen(f"sed -n '1p' {N_raw_stats} | grep -o '^[0-9]*'").read())

        p_N = n / n_total_N
        seed = seed_dict[f"{wildcards.id}_T_reads{wildcards.reads}M"]
        shell("samtools view --no-PG -@ 4 -s " + str(seed + p_N) + " -b -o {output.receptor} {input.N}")
        shell("samtools index {output.receptor}")
        shell("posmix -1 {output.receptor} -2 {input.T} -p {input.posfreq} -s " + str(seed) + " -o {output.posmix_raw}")
        shell('gatk --java-options "-XX:ParallelGCThreads=8" MarkDuplicates -I {output.posmix_raw} -O {output.posmix} -M {output.posmix}.dupmetrics.txt')
        shell("samtools index {output.posmix}")

rule tidy_N:
    input:
        bam = "posmix/{id}_N.bam",
        cp = "cp.txt"
    output:
        reheader = temp("posmix/{id}_N.reheader.bam"),
        tidy = "posmix/{id}_N.tidy.bam"
    threads: 8
    benchmark: "logs/{id}_tidy_N.benchmark.txt"
    shell:
        '''
        samtools view -H {input.bam} | sed '/^@PG/d' | samtools reheader --no-PG /dev/stdin {input.bam} \
        | samtools addreplacerg --no-PG -r '@RG\tID:normal\tSM:normal' --threads 8 -o {output.reheader} /dev/stdin
        
        gatk --java-options "-XX:ParallelGCThreads=8" MarkDuplicates -I {output.reheader} -O {output.tidy} -M {output.tidy}.dupmetrics.txt

        samtools index {output.tidy}
        samtools stats -@ 8 {output.tidy} > {output.tidy}.stat.txt
        '''

rule N_pileup:
    input:
        tidy = "posmix/{id}_N.tidy.bam",
        bed = config["bed_file"]
    output:
        pileup = "varscan/{id}_N.pileup.gz"
    threads: 10
    benchmark: "logs/{id}_N_pileup.benchmark.txt"
    params:
        '-d 1000000 -q 30 -Q 30 -A -B -x --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -f ' + genome_dir
    shell:
        "cat {input.bed} | parallel -k --pipe -j {threads} -N 100 'samtools mpileup {params} --positions /dev/stdin {input.tidy}' | bgzip -c > {output.pileup}"

rule T_pileup:
    input:
        "posmix/{id}_T_reads{reads}M.bam",
        config["bed_file"]
    output:
        "varscan/{id}_T_reads{reads}M.pileup.gz"
    benchmark: "logs/{id}_reads{reads}M_T_pileup.benchmark.txt"
    params:
        '-d 1000000 -q 30 -Q 30 -A -B -x --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -f ' + genome_dir
    threads: 10
    shell:
        '''
        cat {input[1]} | parallel -k --pipe -j {threads} -N 100 'samtools mpileup {params} --positions /dev/stdin {input[0]}' | bgzip -c > {output}
        '''

rule varscan:
    input:
        T = "varscan/{id}_T_reads{reads}M.pileup.gz",
        N = lambda wildcards: N_pileup[wildcards.id]
    output:
        snp = 'varscan/{id}_reads{reads}M.raw.snp.vcf.gz',
        indel = 'varscan/{id}_reads{reads}M.raw.indel.vcf.gz'
    benchmark: "logs/{id}_reads{reads}M_varscan.benchmark.txt" 
    params:
        lcr = os.path.join("/biobase/data/genome", config["genome"], 'lcr.bed'),
        opts = "--min-var-freq " + str(config['min_freq']) + ' --min-coverage-normal 10 --min-coverage-tumor 100 --min-freq-for-hom 0.75 --normal-purity 1.00 --tumor-purity 1.00 --p-value 0.99 --somatic-p-value 0.01 --strand-filter 0 --validation 1 --output-vcf 1'
    shell:
        '''
        opref={output.snp}
        opref=${{opref%.snp.vcf.gz}}
        java -jar VarScan.v2.4.3.jar somatic <(bgzip -d -c {input.N}) <(bgzip -d -c {input.T}) ${{opref}}.all {params.opts}
        cat ${{opref}}.all.snp.vcf | perl /biobase/apps/novotumor/somaticfilter.pl - | bgzip -c > {output.snp}
        cat ${{opref}}.all.indel.vcf | perl /biobase/apps/novotumor/somaticfilter.pl - | bgzip -c > {output.indel}
        bgzip ${{opref}}.all.snp.vcf
        bgzip ${{opref}}.all.indel.vcf
        bgzip ${{opref}}.all.validation
        '''

rule varscan_novo:
    input:
        bam = "posmix/{id}_T_reads{reads}M.bam",
        snp = "varscan/{id}_reads{reads}M.raw.snp.vcf.gz",
        indel = "varscan/{id}_reads{reads}M.raw.indel.vcf.gz"
    output:
        snp_rc = 'varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.gz',
        indel_rc = 'varscan_novo/{id}_reads{reads}M.raw.indel.vcf.rc.gz',
        snp = 'varscan_novo/{id}_reads{reads}M.snp.pass.vcf',
        indel = 'varscan_novo/{id}_reads{reads}M.indel.pass.vcf',
        snp2 = 'varscan_novo/{id}_reads{reads}M.snp.failed.vcf',
        indel2 = 'varscan_novo/{id}_reads{reads}M.indel.failed.vcf'
    benchmark: "logs/{id}_reads{reads}M_varscan_novo.benchmark.txt"
    params:
        fpfilter = "--min-var-freq " + str(config['min_freq']) + ' --min-var-count 4 --min-strandedness 0.05 --min-ref-readpos 0.2 --min-var-readpos 0.2 --min-var-dist3 0.15 --min-ref-dist3 0.2  --min-var-basequal 30 --min-ref-basequal 30 --max-rl-diff 0.2 --max-mmqs-diff 1000000 --max-ref-mmqs 50 --min-ref-mapqual 20 --min-var-mapqual 30',
        ref = genome_dir
    threads: 2
    shell:
        '''
        zcat {input.snp} | awk -F"\\t" '{{if($0!~/^#/){{print $1"\\t"$2"\\t"$2}}}}' | parallel -k -j {threads} --pipe -N 1000 'bam-readcount -q 30 -b 30 -f {params.ref} -l /dev/stdin -w 0 {input.bam}' | bgzip -c > {output.snp_rc}
        zcat {input.indel} | awk -F"\\t" '{{if($0!~/^#/){{if(length($4)>1){{print $1"\\t"$2+1"\\t"$2+1}} else {{print $1"\\t"$2"\\t"$2}}}}}}' | parallel -k -j {threads} --pipe -N 1000 'bam-readcount -q 30 -b 30 -f {params.ref} -l /dev/stdin -w 0 {input.bam}' | bgzip -c > {output.indel_rc}
        java -jar VarScan.v2.4.3.jar fpfilter <(zcat {input.snp}) <(bgzip -d -c {output.snp_rc}) {params.fpfilter} --output-file {output.snp} --filtered-file {output.snp2}
        java -jar VarScan.v2.4.3.jar fpfilter <(zcat {input.indel}) <(bgzip -d -c {output.indel_rc}) {params.fpfilter} --output-file {output.indel} --filtered-file {output.indel2}
        '''

rule bgzip:
    input:
        snp = 'varscan_novo/{id}_reads{reads}M.snp.pass.vcf',
        indel = 'varscan_novo/{id}_reads{reads}M.indel.pass.vcf',
        snp2 = 'varscan_novo/{id}_reads{reads}M.snp.failed.vcf',
        indel2 = 'varscan_novo/{id}_reads{reads}M.indel.failed.vcf'
    output:
        snp_gz = 'varscan_novo/{id}_reads{reads}M.snp.pass.vcf.gz',
        indel_gz = 'varscan_novo/{id}_reads{reads}M.indel.pass.vcf.gz',
        snp2_gz = 'varscan_novo/{id}_reads{reads}M.snp.failed.vcf.gz',
        indel2_gz = 'varscan_novo/{id}_reads{reads}M.indel.failed.vcf.gz'
    benchmark: "logs/{id}_reads{reads}M_bgzip.benchmark.txt"
    shell:
        '''
        bgzip {input.snp}
        bgzip {input.snp2}
        bgzip {input.indel}
        bgzip {input.indel2}
        '''

rule manta:
    input:
        T = "posmix/{id}_T_reads{reads}M.bam",
        N = lambda wildcards: N_tidy[wildcards.id]
    output:
        mantaindel = "manta/{id}_reads{reads}M/results/variants/candidateSmallIndels.vcf.gz",
        vcf = "manta/{id}_reads{reads}M/results/variants/somaticSV.vcf.gz"
    benchmark: "logs/{id}_reads{reads}M_manta.benchmark.txt"
    params:
        outdir = "manta/{id}_reads{reads}M",
        bed = config["bed_file"] + '.gz',
        ref = genome_dir
    threads: 20
    shell:
        '''
        mkdir -p {params.outdir}
        if [ ! -f {params.outdir}/runWorkflow.py ]; then
            configManta.py --outputContig --generateEvidenceBam --exome --normalBam {input.N} --tumorBam {input.T} --referenceFasta {params.ref} --runDir {params.outdir}
        fi
            {params.outdir}/runWorkflow.py -m local -j {threads}
        '''

rule strelka:
    input:
        T = "posmix/{id}_T_reads{reads}M.bam",
        N = lambda wildcards: N_tidy[wildcards.id],
        mantaindel = 'manta/{id}_reads{reads}M/results/variants/candidateSmallIndels.vcf.gz',
        bed_file = config["bed_file"]
    output:
        snv = 'strelka/{id}_reads{reads}M/results/variants/somatic.snvs.vcf.gz',
        indel = 'strelka/{id}_reads{reads}M/results/variants/somatic.indels.vcf.gz'
    benchmark: "logs/{id}_reads{reads}M_strelka.benchmark.txt"
    params:
        outdir= 'strelka/{id}_reads{reads}M',
        ref = genome_dir
    threads: 40
    shell:
        '''
        mkdir -p {params.outdir}
        if [ ! -f {params.outdir}/runWorkflow.py ]; then
            configureStrelkaSomaticWorkflow.py --indelCandidates {input.mantaindel} --callRegions {input.bed_file}.gz --exome --normalBam {input.N} --tumorBam {input.T} --referenceFasta {params.ref} --runDir {params.outdir} --reportEVSFeatures
        fi
        {params.outdir}/runWorkflow.py -m local -j {threads}
        '''

rule bed2intervals:
    input:
        bed = config["bed_file"]
    output:
        outdir = directory("mutect/interval_files/"),
        intervals_list = 'mutect/intervals.list'
    params:
        split = bed_nsplit,
        ref = genome_dir
    shell:
        '''
        mkdir -p {output.outdir}
        cat {input.bed} | cut -f1-3 | grep -v '^#' | perl -lane 'print "$F[0]:$F[1]-$F[2]"' > {output.intervals_list}
        gatk SplitIntervals -R {params.ref} -L {output.intervals_list} -O {output.outdir} --scatter-count {params.split}
        '''

rule mutect:
    input:
        T = "posmix/{id}_T_reads{reads}M.bam",
        N = lambda wildcards: N_tidy[wildcards.id],
        intervals_list = 'mutect/intervals.list'
    output:
        vcf = 'mutect/{id}_reads{reads}M.tmpvcf/{split}.vcf.gz',
        vcfstats = "mutect/{id}_reads{reads}M.tmpvcf/{split}.vcf.gz.stats"
    benchmark: "logs/{id}_reads{reads}M_{split}_mutect.benchmark.txt"
    params:
        interval = "mutect/interval_files/{split}-scattered.interval_list",
        ref = genome_dir
    threads: 4
    shell:
        '''
        gatk Mutect2 -R {params.ref} -I {input.T} -I {input.N} -normal normal -O {output.vcf} -L {params.interval}
        '''

rule merge_mutect:
    input:
        vcf = expand("mutect/{{id}}_reads{{reads}}M.tmpvcf/{split}.vcf.gz", split = splitids),
        vcfstats = [f"mutect/{{id}}_reads{{reads}}M.tmpvcf/{split}.vcf.gz.stats" for split in splitids],
        intervals_list = 'mutect/intervals.list'
    output:
        raw = 'mutect/{id}_reads{reads}M.raw.vcf.gz',
        raw_stats = 'mutect/{id}_reads{reads}M.raw.vcf.gz.stats'
    benchmark: "logs/{id}_reads{reads}M_merge_mutect.benchmark.txt"
    shell:
        '''
        gatk MergeVcfs -O {output.raw} $(echo ' {input.vcf}' | sed 's/ / --INPUT /g')
        gatk MergeMutectStats -O {output.raw_stats} $(echo ' {input.vcfstats}' | sed 's/ / --stats /g')
        '''

rule final_filter:
    input:
        raw = 'mutect/{id}_reads{reads}M.raw.vcf.gz',
        intervals_list = 'mutect/intervals.list'
    output:
        vcf = 'mutect/{id}_reads{reads}M.filter.vcf.gz'
    benchmark: "logs/{id}_reads{reads}M_final_filter.benchmark.txt"
    params:
        ref = genome_dir
    shell:
        '''
        gatk FilterMutectCalls -V {input.raw} -O {output.vcf} -L {input.intervals_list} -R {params.ref} --min-reads-per-strand 1 --min-median-base-quality 20 --normal-p-value-threshold 0.01
        '''

rule tsv:
    input:
        "varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.gz"
    output:
        "varscan_novo/{id}_reads{reads}M.raw.snp.vcf.rc.tsv.gz"
    benchmark: "logs/{id}_reads{reads}M_tsv.benchmark.txt"
    params:
        ref = genome_dir
    shell:
        "reshapeRC -i <(bgzip -d -c {input}) -g {params.ref} | bgzip -c > {output}"

rule sort_index:
    input:
        vcf = "varscan/{id}_reads{reads}M.raw.snp.vcf.gz"
    output:
        sort = "varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz"
    benchmark: "logs/{id}_reads{reads}M_sort_index.benchmark.txt"
    shell:
        '''
        zcat {input.vcf} | awk '{{if (substr($1,1,1)=="#") print $0; else print $0 | "sort -k1,1V -k2,2n"}}' | bgzip -c > {output.sort}
        gatk IndexFeatureFile -I {output.sort}
        '''

rule force_calling:
    input:
        T = "posmix/{id}_T_reads{reads}M.bam",
        N = lambda wildcards: N_tidy[wildcards.id],
        intervals_list = 'mutect/intervals.list',
        sort = "varscan/{id}_reads{reads}M.raw.snp.sort.vcf.gz"
    output:
        vcf = 'mutect_a_vraw/{id}_reads{reads}M.tmpvcf/{split}.vcf.gz',
        vcfstats = "mutect_a_vraw/{id}_reads{reads}M.tmpvcf/{split}.vcf.gz.stats"
    benchmark: "logs/{id}_reads{reads}M_{split}_force_calling.benchmark.txt"
    params:
        ref = genome_dir,
        interval = "mutect/interval_files/{split}-scattered.interval_list"
    threads: 4
    shell:
        '''
        gatk Mutect2 -R {params.ref} -I {input.T} -I {input.N} -normal normal -O {output.vcf} -L {params.interval} --alleles {input.sort} --genotype-filtered-alleles
        '''

rule merge_force_calling:
    input:
        vcf = expand("mutect_a_vraw/{{id}}_reads{{reads}}M.tmpvcf/{split}.vcf.gz", split = splitids),
        vcfstats = [f"mutect_a_vraw/{{id}}_reads{{reads}}M.tmpvcf/{split}.vcf.gz.stats" for split in splitids],
        intervals_list = 'mutect/intervals.list'
    output:
        raw = 'mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz',
        raw_stats = 'mutect_a_vraw/{id}_reads{reads}M.raw.vcf.gz.stats'
    benchmark: "logs/{id}_reads{reads}M_merge_force_calling.benchmark.txt"
    shell:
        '''
        gatk MergeVcfs -O {output.raw} $(echo ' {input.vcf}' | sed 's/ / --INPUT /g')
        gatk MergeMutectStats -O {output.raw_stats} $(echo ' {input.vcfstats}' | sed 's/ / --stats /g')
        '''
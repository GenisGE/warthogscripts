# REFERENCE

REF="/davidData/genis/references/pig/SusScrofa11.1.fasta.gz"

F_OUT = 3852
F_IN_PAIRED = 3 # only for paired-end data
#MIN_MQ = 30

# SOFTWARE
BWA="/home/genis/software/bwa-0.7.17/bwa"
BGZIP="/home/genis/software/htslib-1.9/bgzip"
SAMTOOLS="/home/genis/software/samtools-1.9/samtools"
ADAPTERREMOVAL="/home/genis/software/adapterremoval/build/AdapterRemoval"
NGMERGE="/home/genis/software/NGmerge/NGmerge"
FASTQC="/home/genis/software/FastQC/fastqc"
MULTIQC="/home/genis/.local/bin/multiqc"

OUTMAIN="map"
BWA_THREADS = 10
## Kristian did on origianl
## add fastqc after mapping
## replace adapterremoval with ngmerge
## add flag for paired-end mapping to remove all not properly paired.

## GENIS did
## added fastqc_bam (does fastqc not on mapQ30 filtered bams) and put in mapping multiqc report

configfile:"warts_config.yaml"

         
wildcard_constraints:
    sample = "|".join(config.keys()),
    
rule all:
    input:
        os.path.join(OUTMAIN, "multiqc", "report.html")


rule all_fastqc:
    input:
        expand(os.path.join(OUTMAIN, "fastqc_pretrim", "{sample}"), sample=config.keys())


rule fastqc_pre_trim:
    ## https://stackoverflow.com/a/50882104
    input:
        fq1 = lambda wildcards: config[wildcards.sample][0],
        fq2 = lambda wildcards: config[wildcards.sample][1]
    output:
        directory(os.path.join(OUTMAIN, "fastqc_pretrim", "{sample}"))
    threads: 2
    log:
        os.path.join(OUTMAIN, "fastqc_pretrim", "{sample}.log")
    shell: """ 
    mkdir -p {output}
    {FASTQC} -o {output} --nogroup -f fastq -t {threads} {input} &> {log};
"""

rule multiqc_fastqc:
    input:
        expand(os.path.join(OUTMAIN, "fastqc_pretrim", "{sample}"),sample=config.keys())
    output:
        f=os.path.join(OUTMAIN, "fastqc_pretrim_multiqc", "report.html")
    params:
        dirname = lambda wildcards, output: os.path.dirname(output.f),
        base =  lambda wildcards, output: os.path.basename(output.f).replace(".html", "")
    run:
        dirs = set()
        for f in input:
            dirs.add(f)
        dirs_string = " ".join(dirs)
        shell("{MULTIQC} -f -o {params.dirname} -n {params.base} {dirs_string}")

rule bwa_index_fa:
    input:
        REF
    output:
        REF+".amb"
    shell:
        "{BWA} index {input}"

rule samtools_index_fa:
    input:
        REF
    output:
        REF+".fai"
    shell:
        "{SAMTOOLS} faidx {input}"

# rule trim_fq_adapterremoval:
#     input:
#         fq1 = lambda wildcards: config[wildcards.sample][0],
#         fq2 = lambda wildcards: config[wildcards.sample][1],
#     output:
#         pair1 = temp(os.path.join(OUTMAIN, "trim", "{sample}.pair1.truncated")),
#         pair2 = temp(os.path.join(OUTMAIN, "trim", "{sample}.pair2.truncated")),
#         collapsed = temp(os.path.join(OUTMAIN, "trim", "{sample}.collapsed")),
#         singles = temp(os.path.join(OUTMAIN, "trim", "{sample}.singleton.truncated")),
#         discarded = temp(os.path.join(OUTMAIN, "trim", "{sample}.discarded")),
#         trunc = temp(os.path.join(OUTMAIN, "trim", "{sample}.collapsed.truncated"))
#     params:
#         basedir = lambda wildcards, output: os.path.join(os.path.dirname(output.singles), wildcards.sample)
#     log:
#         settings = os.path.join(OUTMAIN, "trim", "{sample}.settings"),
#     threads: 3
#     shell: """
#     {ADAPTERREMOVAL} --file1 {input.fq1} --file2 {input.fq2} --basename {params.basedir} --collapse --threads {threads}
# """


rule trim_fq:
    input:
        fq1 = lambda wildcards: config[wildcards.sample][0],
        fq2 = lambda wildcards: config[wildcards.sample][1],
    output:
        pair1 = temp(os.path.join(OUTMAIN, "trim", "{sample}_1.fastq")),
        pair2 = temp(os.path.join(OUTMAIN, "trim", "{sample}_2.fastq")),
        collapsed = temp(os.path.join(OUTMAIN, "trim", "{sample}")),
    params:
        basedir = lambda wildcards, output: os.path.join(os.path.dirname(output.pair1), wildcards.sample)
    threads: 10
    log:
        settings = os.path.join(OUTMAIN, "trim", "{sample}.settings"),
    shell: """
    {NGMERGE} -1 {input.fq1} -2 {input.fq2} -o {params.basedir} -m 11 -y -g -n {threads} -u 41 -f {params.basedir}
"""

rule bwa_mem_samtools:
    input:
        p1=rules.trim_fq.output.pair1,
        p2=rules.trim_fq.output.pair2,
        bwaindex = rules.bwa_index_fa.output[0],
        samtoolsindex = rules.samtools_index_fa.output[0],
    output:
        bam = temp(os.path.join(OUTMAIN, "bam", "{sample}_paired.bam")),
    params:
        rg = r"@RG\tID:{sample}-paired\tSM:{sample}" ,
        prefix=lambda wildcards, output: output.bam+".temp",
    log:
        os.path.join(OUTMAIN, "bam", "{sample}.bwa.paired.log"),
    threads: BWA_THREADS
    shell: """
    # -M marks short split hits as secondary
    {BWA} mem -M -R '{params.rg}' -t {threads} {REF} {input.p1}  {input.p2} 2> {log} | {SAMTOOLS} sort -n -T {params.prefix}_n -O BAM -@ 2 | {SAMTOOLS} fixmate -@ 2 -m -O BAM /dev/stdin /dev/stdout | {SAMTOOLS} sort -T {params.prefix}_c -OBAM -@ 2 |  {SAMTOOLS}  calmd -Q - {REF} | {SAMTOOLS} markdup -@ 2 -T  {params.prefix}_markdup /dev/stdin /dev/stdout | {SAMTOOLS} view -f {F_IN_PAIRED} -F {F_OUT} -@ 4 -OBAM -o {output.bam} -
"""

## https://www.biostars.org/p/415831/
rule bwa_mem_samtools_collapsed:
    input:
        collapsed=rules.trim_fq.output.collapsed,
        bwaindex = rules.bwa_index_fa.output[0],
        samtoolsindex = rules.samtools_index_fa.output[0],
    output:
        bam = temp(os.path.join(OUTMAIN, "bam", "{sample}_collapsed.bam")),
    params:
        rg = r"@RG\tID:{sample}-collapsed\tSM:{sample}" ,
        prefix=lambda wildcards, output: output.bam+".temp",
    log:
        os.path.join(OUTMAIN, "bam", "{sample}.bwa.collapsed.log"),
    threads: BWA_THREADS
    shell: """
    # -M marks short split hits as secondary
    {BWA} mem -M -R '{params.rg}' -t {threads} {REF} {input.collapsed} 2> {log} | {SAMTOOLS} sort -T {params.prefix}_c -OBAM -@ 2 |  {SAMTOOLS}  calmd -Q - {REF} | {SAMTOOLS} markdup -@ 2 -T  {params.prefix}_markdup /dev/stdin /dev/stdout | {SAMTOOLS} view -@ 4 -F {F_OUT} -OBAM -o {output.bam} -
"""

rule merge_paired_collapsed:
    input:
        expand(os.path.join(OUTMAIN, "bam", "{{sample}}_{t}.bam"), t=["collapsed", "paired"])
    output:
        bam=os.path.join(OUTMAIN, "bam", "{sample}.bam")
    threads: 4
    shell:
        "{SAMTOOLS} merge -O BAM -@ {threads} {output} {input}"
        
        
rule index_bam:
    input:
        bam = rules.merge_paired_collapsed.output.bam
    output:
        rules.merge_paired_collapsed.output.bam + ".bai"
    shell:
        "{SAMTOOLS} index -@ 2 {input}"

rule flagstat_bam:
    input:
        bam = rules.merge_paired_collapsed.output.bam,
        bai = rules.index_bam.output[0]
    output:
        stats = os.path.join(OUTMAIN, "bam", "{sample}.flagstat"),
    shell:
        "{SAMTOOLS} flagstat {input.bam} > {output}"

rule stats_bam:
    input:
        bam = rules.merge_paired_collapsed.output.bam,
        bai = rules.index_bam.output[0]
    output:
        stats = os.path.join(OUTMAIN, "bam", "{sample}.stats")
    shell:
        "{SAMTOOLS} stats -@ 2 -r {REF} {input.bam} > {output}"

rule idxstats_bam:
    input:
        bam = rules.merge_paired_collapsed.output.bam,
        bai = rules.index_bam.output[0]
    output:
        stats = os.path.join(OUTMAIN, "bam", "{sample}.idxstat")
    shell:
        "{SAMTOOLS} idxstats {input.bam} > {output}"


rule fastqc_bam:
     input:
        bam = rules.merge_paired_collapsed.output.bam
     output:
        directory(os.path.join(OUTMAIN, "fastqc_bam", "{sample}"))
     log:
        os.path.join(OUTMAIN, "fastqc_bam", "{sample}.log")
     threads: 2
     shell: """
        mkdir -p {output}
        {FASTQC} -o {output} --nogroup -f bam -t {threads} {input.bam} &> {log};
"""


rule multiqc_stats:
    input:
        expand(os.path.join(OUTMAIN, "bam", "{sample}.{s}"),sample=config.keys(), s=["flagstat", "idxstat", "stats"]),
        expand(os.path.join(OUTMAIN, "trim", "{sample}.settings"),sample=config.keys()),
        expand(os.path.join(OUTMAIN, "fastqc_bam", "{sample}"),sample=config.keys())
    output:
        f=os.path.join(OUTMAIN, "multiqc", "report.html")
    params:
        dirname = lambda wildcards, output: os.path.dirname(output.f),
        base =  lambda wildcards, output: os.path.basename(output.f).replace(".html", "")
    run:
        dirs = set()
        for f in input:
            dirs.add(os.path.dirname(f))
        dirs_string = " ".join(dirs)

        shell("{MULTIQC} -f -o {params.dirname} -n {params.base} {dirs_string}")

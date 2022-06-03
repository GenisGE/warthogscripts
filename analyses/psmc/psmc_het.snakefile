# snakemake to run psmc. can do it using different filtering in the form of bedfiles as wildcards, and differnt minimum allele support to call hets.
# config needs:
#    samples: dic with pairs of sample:pathtobamfile
#    depths: dic with depth thresholds specific to each sample sample:[mindepth, maxdepth]
#    beds: dic with filtername:path_to_bed_with_filters
#    allele_support: list of allele support values to call heterozygotes
#    plot_params: dic with g:generation_time mu:mutation_rate
#    outmain: path to main output directory

import os

# https://github.com/lh3/psmc
VCFUTILS="/home/genis/software/bcftools/misc/vcfutils.pl"
SAMTOOLS="/home/genis/software/samtools-1.9/samtools"
BCFTOOLS="/home/genis/software/bcftools/bcftools"

PSMC_DIR="/home/genis/software/psmc"
FQ2PSMCFA=os.path.join(PSMC_DIR,"utils/fq2psmcfa")
SPLITFA=os.path.join(PSMC_DIR, "utils/splitfa")
PSMC=os.path.join(PSMC_DIR,"psmc")
#PSMC_PLOT=os.path.join(PSMC_DIR,"utils/psmc_plot.pl")
PSMC_PLOT="scripts/plotPsmc.R"
PSMC_PLOT_BOOT="scripts/plotPsmc_bootstrap.R"


OUTMAIN = config["outmain"]
BEDS = config["beds"]
allele_support = config["allele_support"]

NBOOTS=config["bootstraps"]

wildcard_constraints:
    t = "|".join(allele_support),
    sample = "|".join(config["samples"].keys()),
    bed = "|".join(BEDS.keys()),


subworkflow geno_calls:
    workdir:
        "../geno_calls/"
    snakefile:
        "../geno_calls/call.snakefile"
    configfile:
        "../geno_calls/config.yaml"


    
rule all:
    input:
        os.path.join(OUTMAIN, "psmc_output_figures", "all_filters_3_bootstrap.png"),
        expand(os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}.psmc"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),
        expand(os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.png"),
               bed=BEDS.keys(),
               t=allele_support
        ),
        expand(os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}.het"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),
        expand(os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "bootstrap", "{sample}_{i}.psmc"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support,
               i=list(range(1, NBOOTS+1)))





# rules to run psmc

rule gen_fq_mindepthx:
    input:
        
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.fq.gz"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0], # 30/3
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],  # 30*2
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        """



rule gen_psmcfa:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.fq.gz"),
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    shell: """
        {FQ2PSMCFA} -q 20 {input} > {output}
"""



rule split_psmcfa:
    input:
        psmcfa = os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    output:
        split_psmcfa = os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "bootstrap", "{sample}_split.psmcfa"),
    shell:"""
    {SPLITFA} {input.psmcfa} > {output.split_psmcfa}
    """



rule run_psmc:
    input:
        psmcfa = os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    output:
        psmc = os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}.psmc")
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output.psmc} {input.psmcfa}"""



rule bootstrap_psmc:
    input:
        split_psmcfa = os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "bootstrap", "{sample}_split.psmcfa"),
    output:
        psmc = os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "bootstrap", "{sample}_{i}.psmc")
    shell:"""
    {PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -b -o {output.psmc} {input.split_psmcfa}
"""



rule plot_psmc:
    input:
        expand(os.path.join(OUTMAIN, "psmc_output", "{{bed}}", "{{t}}", "{sample}.psmc"),
               sample=config["samples"].keys())
    output:
        png=os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.png")
    params:
        indir = os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}"),
        g=config["plot_params"]["g"],
	mu=config["plot_params"]["mu"]
    shell: "Rscript {PSMC_PLOT} -d {params.indir} -g {params.g} -m {params.mu} -o {output.png}"



rule plot_psmc_bootstrap:
    input:
        expand(os.path.join(OUTMAIN, "psmc_output", "{{bed}}", "{{t}}", "{sample}.psmc"),
               sample=config["samples"].keys()),
        expand(os.path.join(OUTMAIN, "psmc_output", "{{bed}}", "{{t}}", "bootstrap", "{sample}_{i}.psmc"),
               sample=config["samples"].keys(), i  = list(range(1, NBOOTS+1)))
    output:
        png=os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}_bootstrap.png")
    params:
        indir = os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}"),
        bootdir = os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "bootstrap"),
        g=config["plot_params"]["g"],
	mu=config["plot_params"]["mu"]
    shell: "Rscript {PSMC_PLOT_BOOT} -d {params.indir} -b {params.bootdir} -g {params.g} -m {params.mu} -o {output.png}"

 


# rules to estimate genome_wide heterozygoisties with bcftools, using same sites as in psmc

rule get_het_mindepthx:
    input:
        geno_calls(os.path.join("results", "vcf", "{sample}.bcf.gz"))
    output:
        os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}.bcf.stats"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0], # 30/3
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],  # 30*2
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk | {BCFTOOLS} stats -s - > {output}
        """


rule estimate_het:
    input:
       os.path.join(OUTMAIN, "bcf_stats", "{bed}", "{t}", "{sample}.bcf.stats")
    output:
       os.path.join(OUTMAIN, "hets", "{sample}_{bed}_{t}.het")
    shell: "grep '^PSC' {input} | awk '{{print $6/($4+$5+$6)}}' > {output}"



BCFTOOLS="/newHome/genis/mydata/software/bcftools/bcftools"
BGZIP="/newHome/genis/mydata/software/htslib-1.9/bgzip"
TABIX="/newHome/genis/mydata/software/htslib-1.9/tabix"
RSCRIPT="Rscript"

SCRIPTSDIR=config["scriptsdir"]
TT=os.path.join(SCRIPTSDIR, "tt.R")
TTO=os.path.join(SCRIPTSDIR, "tto.R")
PLOT=os.path.join(SCRIPTSDIR, "collectPlotPointEstimates.R")

OUTMAIN=config["outmain"]



rule all_tt:
    input:
        expand(os.path.join(OUTMAIN, "estimates", "tt", "{pair}.tt..{f}"),
               pair=config["do_tt"].keys(),
               f=["params.res", "split.years"]
        ),


rule all_tto:
    input:
        expand(os.path.join(OUTMAIN, "estimates", "tto", "{pair}.tto.{f}"),
               pair= config["do_tto"].keys(),
               f=["params.res", "split.years"]
        )


rule select_polarized_sites:
    input:
        bcf = config["bcf"]
    output:
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz")
    params:
        outgroupsamples = ",".join(config["outgroups"])
    shell:
        """{BCFTOOLS} view -s {params.outgroupsamples} {input.bcf} | {BCFTOOLS} view  -e 'GT=="1/1" || GT=="0/1" || GT=="./."'  | {BCFTOOLS} query -f '%CHROM\t%POS\n' | {BGZIP} -c > {output.sites}"""


rule index_polarized_sites:
    input:
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz")
    output:
        index = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz.tbi")
    shell:
        """{TABIX} -b 2 -e 2 {input.sites}"""


rule select_tto_sites:
    input:
        bcf = config["bcf"],
        polar_sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz")
    output:
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites_outgroupderived.txt.bgz")
    params:
        close_outgroup = config["close_outgroup"]
    shell:
        """{BCFTOOLS} view -s {params.close_outgroup} {input.bcf} | {BCFTOOLS} view  -i 'GT=="1/1"' |  {BCFTOOLS} query -f '%CHROM\t%POS\n' | {BGZIP} -c > {output.sites}"""



rule index_tto_sites:
    input:
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites_outgroupderived.txt.bgz")
    output:
        index = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites_outgroupderived.txt.bgz.tbi")
    shell:
        """{TABIX} -b 2 -e 2 {input.sites}"""



rule get_genomewide_2dsfs:
    input:
        bcf = config["bcf"],
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz"),
        sites_index = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites.txt.bgz.tbi"),
    params:
        s1 = lambda wildcards: config["do_tt"][wildcards.pair][0],
        s2 = lambda wildcards: config["do_tt"][wildcards.pair][1],
    output:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pair}_unfolded.2dsfs")
    shell:
        """{BCFTOOLS} view -s {params.s1},{params.s2} -T {input.sites} {input.bcf} | {BCFTOOLS} view  -e 'GT=="./."'  | {BCFTOOLS} query -f "[%GT ]\n" | awk '{{a[$1"_"$2]++}} END{{for(g in a){{print g, a[g]}}}}' > {output.sfs}"""



rule get_genomewide_2dsfs_outgroupderived:
    input:
        bcf = config["bcf"],
        sites = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites_outgroupderived.txt.bgz"),
        sites_index = os.path.join(OUTMAIN, "sites", "ancestral_supported_sites_outgroupderived.txt.bgz.tbi"),
    params:
        s1 = lambda wildcards: config["do_tt"][wildcards.pair][0],
        s2 = lambda wildcards: config["do_tt"][wildcards.pair][1],
    output:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pair}_unfolded_outgroupderived.2dsfs")
    shell:
        """{BCFTOOLS} view -s {params.s1},{params.s2} -T {input.sites} {input.bcf} | {BCFTOOLS} view  -e 'GT=="./."'  | {BCFTOOLS} query -f "[%GT ]\n" | awk '{{a[$1"_"$2]++}} END{{for(g in a){{print g, a[g]}}}}' > {output.sfs}"""




rule get_tt_estimates:
    input:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pair}_unfolded.2dsfs")
    output:
        os.path.join(OUTMAIN, "estimates", "tt", "{pair}.tt.params.res"),
        os.path.join(OUTMAIN, "estimates", "tt", "{pair}.tt.split.years"),
    params:
        outprefix = os.path.join(OUTMAIN, "estimates", "tt", "{pair}"),
        g = config["g"],
        mu = config["mu"]
    shell:
        """
        {RSCRIPT} {TT} {input.sfs} {params.outprefix} {params.mu} {params.g}
"""



rule get_tto_estimates:
    input:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{pair}_unfolded_outgroupderived.2dsfs")
    output:
        os.path.join(OUTMAIN, "estimates", "tto", "{pair}.tto.params.res"),
        os.path.join(OUTMAIN, "estimates", "tto", "{pair}.tto.split.years"),
    params:
        outprefix = os.path.join(OUTMAIN, "estimates", "tto", "{pair}"),
        g = config["g"],
        mu = config["mu"]
    shell:
        """
        {RSCRIPT} {TTO} {input.sfs} {params.outprefix} {params.mu} {params.g}
"""




rule plot_point_estimates:
    input:
        lambda wildcards: expand(os.path.join(OUTMAIN, "estimates", "{{method}}", "{pair}.{{method}}.{f}"),
                                 pair=config["do_{m}".format(m=wildcards.method)].keys(),
                                 f=["params.res", "split.years"])
    output:
        png = os.path.join(OUTMAIN, "plots", "split_times_{method}.png")
    params:
        indir = os.path.join(OUTMAIN, "estimates", "{method}")
    shell:
        """
        {RSCRIPT} {PLOT} {params.indir} {output.png} {wildcards.method} 
"""


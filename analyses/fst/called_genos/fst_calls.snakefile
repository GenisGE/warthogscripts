import itertools as it
import pandas as pd


BCFTOOLS="/home/genis/software/bcftools/bcftools"
#ANGSD="/home/genis/software/angsd"
#RM_INDELS_AWK="/home/genis/impala/analyses/goatMapV2/psmc/scripts/rm_indels.awk" 
R="Rscript"
REICHFST="scripts/reich_fst.R"


OUTMAIN=config["outmain"]
BEDS = config["beds"]
OUTMAIN = config["outmain"]
BEDS = config["beds"]
allele_support = config["allele_support"]
AUTOSOMES=config["autosomes"]

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
        os.path.join(OUTMAIN,"collected_fst.txt"),
        os.path.join(OUTMAIN, "formatted_fst.tsv") # the rule that produces this requries a popinfo.list file where first column are population and second column sample id



rule clean_bcf:
    input:
       bcf = geno_calls(os.path.join("results", "vcf", "{sample}.bcf.gz")),
       idx = geno_calls(os.path.join("results", "vcf", "{sample}.bcf.gz.csi"))
    output:
       bcf = temp(os.path.join(OUTMAIN, "input", "{sample}_{t}_{bed}.bcf.gz")),
       idx = temp(os.path.join(OUTMAIN, "input", "{sample}_{t}_{bed}.bcf.gz.csi"))
    params:
       bed = lambda wildcards: BEDS[wildcards.bed],
       mindepth = 10
    shell: """
      {BCFTOOLS} view -R {AUTOSOMES} -T {params.bed} -i 'sum(INFO/DP4)>={params.mindepth}' -V mnps,indels -Ou {input.bcf} | {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t}) || GT=="hom"' -Ob -o {output.bcf}
      {BCFTOOLS} index {output.bcf}
"""

    

rule sfs2d_bcftools:
    input:
        bcf1 = os.path.join(OUTMAIN, "input", "{s1}_{t}_{bed}.bcf.gz"),
        bcf2 = os.path.join(OUTMAIN, "input", "{s2}_{t}_{bed}.bcf.gz"),
	idx1 = os.path.join(OUTMAIN, "input", "{s1}_{t}_{bed}.bcf.gz.csi"),
	idx2 = os.path.join(OUTMAIN, "input", "{s2}_{t}_{bed}.bcf.gz.csi")
    output:
        os.path.join(OUTMAIN, "sfs_2d", "{s1}_{s2}_bcftools_{t}_{bed}.sfs")
    params:
        regions = AUTOSOMES
    shell: """
       {BCFTOOLS} merge -R {params.regions} {input.bcf1} {input.bcf2} | bcftools query -f "[%GT ]\n" | awk '{{a[$1"_"$2]++}} END{{for(g in a){{print g, a[g]}}}}' > {output}
"""




rule collect_2dsfs_bcftools:
    input:
       expand(os.path.join(OUTMAIN, "sfs_2d", "{s[0]}_{s[1]}_bcftools_{t}_{bed}.sfs"), s=it.combinations(config["samples"].keys(), 2), t=allele_support, bed=BEDS.keys())
    output:
       f = os.path.join(OUTMAIN,"collected_2dsfs.txt")
    run:
        import os
        import pandas as pd
	import re
	keys = ['0/0_0/0', '0/0_0/1', '0/0_1/1', '0/1_0/0', '0/1_0/1', '0/1_1/1', '1/1_0/0', '1/1_0/1', '1/1_1/1']
        data = []
        names = []
        for x in input:
            #name = re.sub("_bcftools.*.sfs", "", os.path.basename(x))
            name = re.sub(".sfs", "", os.path.basename(x))
            names.append(name)
        a = pd.DataFrame(index=names, columns=["aaAA","aaAD","aaDD","adAA","adAD","adDD","ddAA","ddAD","ddDD"])
        for x in input:
            d={}
            with open(x, 'r') as fh:
                name = re.sub(".sfs", "", os.path.basename(x))
                for t in fh.readlines():
                    l = t.strip().split(" ")
                    d[l[0]] = int(l[1])
            a.loc[name] = [d[k] for k in keys]
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")



rule get_fst_reich:
    input:
        f=os.path.join(OUTMAIN, "collected_2dsfs.txt")
    output:
        f = os.path.join(OUTMAIN, "collected_fst.txt")
    shell:"""
    {R} {REICHFST} {input.f} {output.f}
"""


rule format_fst_tables:
    input:
        f = os.path.join(OUTMAIN, "collected_fst.txt")
    output:
        f = os.path.join(OUTMAIN, "formatted_fst.tsv")
    params:
        popinfo = "popinfo.list"
    run:
        popinfo = pd.read_table(params.popinfo, header=None)
        popinfo.index = popinfo[1].apply(str)
        fst = pd.read_table(input.f, header=None, sep=" ")
        newdf = pd.DataFrame()
        
        newdf["id1"] = fst[0].apply(lambda x: x.split("_")[0])
        newdf["id2"] = fst[0].apply(lambda x: x.split("_")[1])
        newdf["pop1"] = popinfo[0].loc[newdf["id1"]].values
        newdf["pop2"] = popinfo[0].loc[newdf["id2"]].values
        newdf["ad"] = fst[0].apply(lambda x: x.split("_")[3])
        newdf["fst"] = fst[1]

        newdf.to_csv(output.f, sep="\t", header=True, index=False)





        

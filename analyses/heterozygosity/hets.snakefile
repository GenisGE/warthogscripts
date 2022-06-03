# stolen and modified from kristian snakemake in /home/leopard/users/krishang/old/relatedness/run_2dsfs_cleanref.snakefile
# this only only estiamtes heterozygosity from sfs (not pairwise ibs relatedness)
# config needs:
#   ref: path to reference gneome. Will be used as ancestral to polarize saf (doesnt matter for heterozygositie if is not really ancestral)
#   outmain: path to main output folder
#   angsdparams: dict with  minQ: minimum baseq uality minMapQ: minimum mappingquality
#   chroms: paht to file with chromosomes to use (-rf)
#   sites: path tho agnsd format regions file with regions to use (-sites)
#   samples: dict where key:value are sample_id:pathtobam

import itertools as it

ANGSD="/home/genis/software/angsd"

OUTMAIN=config["outmain"]
SAMPLES = [x for x in config["samples"].keys()]


rule all:
    input:
        os.path.join(OUTMAIN, "sfs", "collected.txt"),

        
rule per_sample_saf:
    input:
        bam= lambda wildcards: config["samples"][wildcards.s]
    output:
        saf_idx = os.path.join(OUTMAIN, "safs", "{s}.saf.idx"),
        saf = os.path.join(OUTMAIN, "safs", "{s}.saf.gz"),
        saf_pos = os.path.join(OUTMAIN, "safs", "{s}.saf.pos.gz"),
        arg = os.path.join(OUTMAIN, "safs", "{s}.arg"),
    params:
        outbase = lambda wildcards, output: output.saf_idx.replace(".saf.idx", ""),
        minQ = config["angsdparams"]["minQ"],
        minMapQ = config["angsdparams"]["minMapQ"],
        chroms = config["chroms"],
        sites = config["sites"],
        anc = config["ref"]
    threads: 3
    shell:
        "{ANGSD}/angsd -i {input.bam} -out {params.outbase} -minQ {params.minQ} -minMapQ {params.minMapQ} -dosaf 1 -sites {params.sites} -rf {params.chroms} -anc {params.anc} -GL 2"

        
rule sfs:
    input:
        saf_idx1 = os.path.join(OUTMAIN, "safs", "{s}.saf.idx"),
    output:
        os.path.join(OUTMAIN, "sfs", "{s}.sfs")
    threads: 15
    log:
        os.path.join(OUTMAIN, "sfs", "{s}.pre_log")
    shell:
        "{ANGSD}/misc/realSFS -P {threads} {input.saf_idx1} > {output} 2> {log}"



rule clean_sfs_log:
    input:
        log = os.path.join(OUTMAIN, "sfs", "{s}.pre_log")
    output:
        log = os.path.join(OUTMAIN, "sfs", "{s}.log")
    shell: "tail -100 {input.log} > {output.log}"



rule collect_sfs:
    input:
        sfs = expand(os.path.join(OUTMAIN, "sfs", "{s}.sfs"), s=SAMPLES),
        clean_logs = expand(os.path.join(OUTMAIN, "sfs", "{s}.log"), s=SAMPLES)
    output:
        f=os.path.join(OUTMAIN, "sfs", "collected.txt")
    run:
        import os
        import pandas as pd
        data = []
        names = []
        for x in input.sfs:
            name = os.path.basename(x).replace(".sfs", "")
            names.append(name)
            with open(x, 'r') as fh:
                t = fh.read()
                data.append([float(x) for x in t.rstrip().split()])
        a = pd.DataFrame(data, index=names, columns=["aa","ad","dd"])
        a["het"] = a["ad"] / a.sum(1)
        a.to_csv(output.f, index=True, header=True, index_label="id", sep=" ")



# snakefile to do abbababa1. This one uses anc as outgroup, and uses all bases (no snp calling).

import pandas as pd
import re

ANGSD="/home/genis/software/angsd/angsd"
JACKKNIFE="/home/genis/software/angsd/R/jackKnife.R"

OUTMAIN="results"

chromfile=config["chroms"]
with open(chromfile, "r") as fh:
    CHROMS=[x.rstrip() for x in fh.readlines()]


wildcard_constraints:
    chrr = "|".join(CHROMS)



rule all:
    input:
        os.path.join(OUTMAIN, "dstats", "out_dstats.txt")


rule do_idlist:
    input:
        bamlist = config['bamlist']
    output:
        ids =  os.path.join(OUTMAIN, "indlists", "list.ids")
    run:
        with open(input.bamlist, "r") as fh:
            bams = [x.rstrip() for x in fh.readlines()]
            #bams = bams[:-1] this only with -useLast 1
        with open(output.ids, "w+") as fh:
            for b in bams:
                s = os.path.basename(b).replace(".bam", "")
                fh.write(f'{s}\n')

                
rule do_abbababa_per_chrom:
    input:
        bamlist = config["bamlist"]
    output:
        temp(os.path.join(OUTMAIN, "abbababa", "{chrr}.abbababa"))
    params:
        outprefix = os.path.join(OUTMAIN, "abbababa", "{chrr}"),
        r = "{chrr}",
        regions = config['bed'],
        MINQ = 30,
        MINMAPQ = 30,
        blocksize = 5000000,
        anc = config['anc']
    log: os.path.join(OUTMAIN, "abbababa", "{chrr}.arg")
    shell: "{ANGSD} -doAbbababa 1 -doCounts 1 -GL 2 -out {params.outprefix}  -bam {input.bamlist} -minmapQ {params.MINMAPQ} -minq {params.MINQ} -r {params.r} -sites {params.regions} -blockSize {params.blocksize} -anc {params.anc} -useLast 0"




rule concat_abbababa:
    input:
        expand(os.path.join(OUTMAIN, "abbababa", "{chrr}.abbababa"), chrr=CHROMS)
    output:
        f = os.path.join(OUTMAIN, "abbababa", "out.abbababa")
    run:
        f1 = input[0]
        shell("head -1 {f1} > {output.f}")
        for f in input:
            shell("cat {f} | sed 1d >> {output.f}")



rule do_jackknife:
    input:
        abbababa = os.path.join(OUTMAIN, "abbababa", "out.abbababa"),
        ids = os.path.join(OUTMAIN, "indlists", "list.ids")
    output:
        os.path.join(OUTMAIN, "dstats", "out_dstats.txt")
    params:
        outprefix = os.path.join(OUTMAIN, "dstats", "out_dstats")
    shell: "Rscript {JACKKNIFE} file={input.abbababa} indNames={input.ids} outfile={params.outprefix}"
        

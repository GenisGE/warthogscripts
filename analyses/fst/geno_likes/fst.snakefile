# snakemake to estimate 2dsfs, sfs and fst between popualtions pairs.
# version 2 intention is to do a more efficient version which will use downsampled saf
# reminder to explain what config needs. also put angsd parameters to config would be nice


import pandas as pd
import itertools as it

ANGSD="/home/genis/software/angsd/angsd"
REALSFS="/home/genis/software/angsd/misc/realSFS"

#OUTMAIN="/emc/genis/impala/sfss_fst/impalaMap/results"

OUTMAIN=config["outmain"]
OUTBIG=config["outbig"]

infofile=config["info"]


info = pd.read_table(infofile, header=None, names=["pops", "bams"])
pops = list(info["pops"].unique())
pops.remove("Zambia") # hardcoded for warthogs

rule all_saf:
    input:
        expand(os.path.join(OUTMAIN, "safs", "{p}.saf.gz"), p=pops)
    
rule all:
    input:
        expand(os.path.join(OUTMAIN, "sfs", "{p}.sfs"), p=pops),
        os.path.join(OUTMAIN, "fst", "all_fst.txt"),
#        expand(os.path.join(OUTMAIN, "pbs", "{trio}_slidingwindow_{win}_stepsize_{step}.pbs"), trio=config["trios"].keys(),
#               win=config["windowparams"]["win"], step = config["windowparams"]["step"])
#    shell: "ln -s -f {OUTMAIN} ."


        
rule do_bam_lists:
    input:
        ancient(config['info'])
    output:
        expand(os.path.join(OUTMAIN, "bamlists", "bams_{p}.list"), p = pops),
        #expand(os.path.join(OUTMAIN, "bamlists", "samples_{pop}.list"), pop = pops)
    run:
        for p in pops:
            bams = info.loc[info["pops"]==p].bams
            outfile1 = os.path.join(OUTMAIN, "bamlists", "bams_{}.list".format(p))
            bams.to_csv(outfile1, header=False, index = False)

            #ids = info.loc[info.my_locality2==pop].Serial_number
            #outfile2 = os.path.join(OUTMAIN, "bamlists", "samples_{}.list".format(pop))
            #ids.to_csv(outfile2, header=False, index = False)



rule do_saf_full:
    input:
        bamlist = ancient(os.path.join(OUTMAIN, "bamlists", "bams_{p}.list"))
    output:
        saf = os.path.join(OUTBIG, "safs", "{p}.saf.gz"),
        saf_idx = os.path.join(OUTBIG, "safs", "{p}.saf.idx"),
        saf_pos = os.path.join(OUTBIG, "safs", "{p}.saf.pos.gz")
    params:
        anc = config['anc'],
        outprefix = os.path.join(OUTBIG, "safs", "{p}"),
        sites = config['sites'],
        rf = config['chroms'],
        minQ = 30,
        minMapQ = 30
    log: os.path.join(OUTBIG, "safs", "{p}.arg")
    threads: 5
    shell: "{ANGSD} -b {input.bamlist} -gl 2 -dosaf 1 -anc {params.anc} -out {params.outprefix} -P {threads} -rf {params.rf} -sites {params.sites} -minQ {params.minQ} -minMapQ {params.minMapQ} 2> /dev/null"

           
rule link_saf:
    input:
        saf = ancient(rules.do_saf_full.output.saf),
        saf_idx = ancient(rules.do_saf_full.output.saf_idx),
        saf_pos = ancient(rules.do_saf_full.output.saf_pos)
    output:
        saf = os.path.join(OUTMAIN, "safs", "{p}.saf.gz"),
        saf_idx = os.path.join(OUTMAIN, "safs", "{p}.saf.idx"),
        saf_pos = os.path.join(OUTMAIN, "safs", "{p}.saf.pos.gz")
    shell: """
    ln -s -f {input.saf} {output.saf}
    ln -s -f {input.saf_idx} {output.saf_idx}
    ln -s -f {input.saf_pos} {output.saf_pos}
"""

           
rule do_saf_small:
    input:
        bamlist = os.path.join(OUTMAIN, "bamlists", "bams_{p}.list")
    output:
        saf = os.path.join(OUTMAIN, "small_safs", "{p}.saf.gz"),
        saf_idx = os.path.join(OUTMAIN, "small_safs", "{p}.saf.idx"),
        saf_pos = os.path.join(OUTMAIN, "small_safs", "{p}.saf.pos.gz")
    log: os.path.join(OUTMAIN, "small_safs", "{p}.arg")         
    params:
        sites = config['small_sites'],
        rf = config['small_chroms'],
        anc = config['anc'],
        outprefix = os.path.join(OUTMAIN, "small_safs", "{p}"),
        minQ = 30,
        minMapQ = 30
    threads: 5
    shell: "{ANGSD} -b {input.bamlist} -gl 2 -dosaf 1 -anc {params.anc} -out {params.outprefix} -P {threads} -rf {params.rf} -sites {params.sites} -minQ {params.minQ} -minMapQ {params.minMapQ} 2> /dev/null"
        
    
rule do_sfs:
    input:
        idx = os.path.join(OUTMAIN, "small_safs", "{p}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "sfs", "{p}.sfs")
    threads: 15
    log: os.path.join(OUTMAIN, "sfs", "{p}.log")
    shell: "{REALSFS} {input.idx} -P {threads} > {output.sfs} 2> {log}"



rule do_2dsfs:
    input:
        idx1 = os.path.join(OUTMAIN, "small_safs", "{p1}.saf.idx"),
        idx2 = os.path.join(OUTMAIN, "small_safs", "{p2}.saf.idx")
    output:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}.sfs")
    log: os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}.log")
    threads: 15
    shell: "{REALSFS} {input.idx1} {input.idx2} -P {threads} > {output.sfs} 2> {log}"


           
rule do_fst_idx:
    input:
        sfs = os.path.join(OUTMAIN, "2dsfs", "{p1}_{p2}.sfs"),
        idx1 = os.path.join(OUTMAIN, "safs", "{p1}.saf.idx"),
        idx2 = os.path.join(OUTMAIN, "safs", "{p2}.saf.idx")
    output:
        fst_idx = temp(os.path.join(OUTMAIN, "fst", "{p1}_{p2}.fst.idx")),
        fst_gz = temp(os.path.join(OUTMAIN, "fst", "{p1}_{p2}.fst.gz"))
    params:
        outprefix = os.path.join(OUTMAIN, "fst","{p1}_{p2}"),
        whichfst = 1
    threads: 15
    shell: "{REALSFS} fst index {input.idx1} {input.idx2} -sfs {input.sfs} -fstout {params.outprefix} -whichFst {params.whichfst} -P {threads} 2> /dev/null"



rule do_fst:
    input:
        fst_idx = os.path.join(OUTMAIN, "fst", "{p1}_{p2}.fst.idx"),
        fst_gz = os.path.join(OUTMAIN, "fst", "{p1}_{p2}.fst.gz")
    output:
        fst = os.path.join(OUTMAIN, "fst", "{p1}_{p2}.fst")
    shell: "{REALSFS} fst stats {input.fst_idx} > {output.fst}"

           

rule collect_fst:
    input:
        ancient(expand(os.path.join(OUTMAIN, "fst", "{p[0]}_{p[1]}.fst"), p = it.combinations(pops, 2)))
    output:
        f = os.path.join(OUTMAIN, "fst", "all_fst.txt")
    run:
        data = []
        for x in input:
            pops = os.path.basename(x).replace(".sfs", "").split("_")
#            popss.append(pops)
            with open(x, 'r') as fh:
                ff = fh.read()
                fst = [float(x) for x in ff.rstrip().split()]
            data.append(pops+fst)
        out = pd.DataFrame(data, columns = ['pop1', 'pop2', 'fst_unweight', 'fst_weight'])
        out.to_csv(output.f, index=False, header=True, sep = "\t")


# rules to run pbs:
         
rule do_fst_index_3pop:
    input:
        sfss = lambda wildcards: expand(os.path.join(OUTMAIN, "2dsfs", "{p[0]}_{p[1]}.sfs"), p=it.combinations(config['trios'][wildcards.trio],2)),
        saf_idxs = lambda wildcards: expand(os.path.join(OUTMAIN, "safs", "{p}.saf.idx"), p=config['trios'][wildcards.trio])
    output:
        fst_idx = os.path.join(OUTBIG, "pbs", "{trio}.fst.idx"),
        fst_gz = os.path.join(OUTBIG, "pbs", "{trio}.fst.gz")
    params:
        whichfst = 1,
        outbase = os.path.join(OUTBIG, "pbs", "{trio}")
    threads: 15
    shell: "{REALSFS} fst index {input.saf_idxs[0]} {input.saf_idxs[1]} {input.saf_idxs[2]} -sfs {input.sfss[0]} -sfs {input.sfss[1]} -sfs {input.sfss[2]} -fstout {params.outbase} -P {threads} 2> /dev/null"

           
rule link_fst_index:
    input:
        fst_idx = ancient(rules.do_fst_index_3pop.output.fst_idx),
        fst_gz = ancient(rules.do_fst_index_3pop.output.fst_gz)
    output:
        fst_idx = os.path.join(OUTMAIN, "pbs", "{trio}.fst.idx"),
        fst_gz = os.path.join(OUTMAIN, "pbs", "{trio}.fst.gz")
    shell: """
    ln -s -f {input.fst_idx} {output.fst_idx}
    ln -s -f {input.fst_gz} {output.fst_gz}
"""

        
rule do_pbs:
    input:
        fst_idx = rules.link_fst_index.output.fst_idx,
    output:
        pbs = os.path.join(OUTMAIN, "pbs", "{trio}_slidingwindow_{win}_stepsize_{step}.pbs")
    params:
        win = "{win}",
        step = "{step}"
    shell: "{REALSFS} fst stats2 {input.fst_idx} -win {params.win} -step {params.step} > {output}"

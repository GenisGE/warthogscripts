# SNAKEMAKE to estimate error rates with 'perfect sample' approach from ANGSD http://www.popgen.dk/angsd/index.php/Error_estimation#Error_estimation_using_an_outgroup_and_an_error_free_individual

# configfile need to have:
# bamlist: list of bam paths from which to estimate error rates
# perfect: path to bam with sample assumed to be error free
# ref: path to fasta reference
# chroms: list chromosomes to use (only chromosome name as in bam file, one per line)
# sites: list of regions to use (with angsd sites format, see http://www.popgen.dk/angsd/index.php/Sites)


SAMTOOLS="/home/genis/software/samtools-1.9/samtools"
ANGSD="/home/genis/software/angsd/angsd"
EST_ERROR="/home/genis/software/angsd/R/estError.R"


OUTMAIN="out_allFilters"


configfile:"config_errrates.yaml"

           
rule all:
    input:
         os.path.join(OUTMAIN, "angsdErrorEst.txt")

rule index_ref:
    input:
        ref = config['ref']
    output:
        idx = config['ref'] + ".fai"
    shell: "{SAMTOOLS} faidx {input.ref}"


rule do_perfect_fasta:
    input:
        bam = config['perfect']
    output:
        fasta = os.path.join(OUTMAIN, os.path.basename(config['perfect']).split('.')[0] + "_errorfree.fa.gz")
    params:
        rf=config['chroms'],
        sites=config['sites'],
        outprefix=os.path.join(OUTMAIN, os.path.basename(config['perfect']).split('.')[0] + "_errorfree")
    log: os.path.join(OUTMAIN, os.path.basename(config['perfect']).split('.')[0] + "_errorfree.arg")
    shell: "{ANGSD} -i {input.bam} -doCounts 1 -doFasta 2 -rf {params.rf} -sites {params.sites} -out {params.outprefix} -minMapQ 30 -minQ 30"

           
rule do_error_rates:
    input:
        idx = rules.index_ref.output.idx,
        perfect = rules.do_perfect_fasta.output.fasta,
        anc=config['ref'],
        bams = config['bamlist']
    output:
        ancerror = os.path.join(OUTMAIN, 'angsderrates.ancError')
    params:
        rf=config['chroms'],
        sites=config['sites'],
        outprefix=os.path.join(OUTMAIN, 'angsderrates')
    log: os.path.join(OUTMAIN, 'angsderrates.arg')
    threads: 10
    shell: "{ANGSD} -doAncError 1 -bam {input.bams} -P {threads} -rf {params.rf} -sites {params.sites} -out {params.outprefix} -minMapQ 30 -minQ 30 -anc {input.anc} -ref {input.perfect}"

        
rule est_error_rates:
    input:
        ancerror = rules.do_error_rates.output.ancerror
    output:
        os.path.join(OUTMAIN, "angsdErrorEst.txt")
    params:
        outprefix = os.path.join(OUTMAIN, "angsdErrorEst")
    shell: "Rscript {EST_ERROR} file={input.ancerror} out={params.outprefix}"


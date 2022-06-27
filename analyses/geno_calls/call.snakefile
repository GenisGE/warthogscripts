##### see http://localhost:8888/notebooks/leopard/misc.ipynb#gen-config-joint_vcf-against-cat

import pandas as pd

BCFTOOLS="/home/genis/software/bcftools/bcftools"
VCF2EIGENSTRAT="python /home/genis/software/gdc/vcf2eigenstrat.py"
#CONVERTER="python3 /home/genis/software/PPP/pgpipe/vcf_format_conversions.py"
PLINK="/home/genis/software/plink"

MIN_MQ=30
MIN_BQ=30
OUTBIG=config["outbig"]
OUTMAIN=config["outmain"]
REF=config["ref"]

## this is a filter used together with total depth and allele depth
#BED = config["bed"]

wildcard_constraints:
    sample = "|".join(config["samples"].keys()),
    chrom="|".join(config["chroms"])


rule all:
    input:
        #multiext(os.path.join(OUTMAIN, "vcf", "merged_snps_filtered"),".bcf.gz", ".bcf.gz.csi"),
        #multiext(os.path.join(OUTMAIN, "eigenstrat", "merged_snps_filtered"), ".ind", ".snp", ".geno"),
        #multiext(os.path.join(OUTMAIN, "plink", "merged_snps_filtered"), ".fam", ".bim", ".bed")



rule do_samples_file:
    output:
        samples_file = os.path.join(OUTBIG, "vcf", "samples.list")
    run:
        s = pd.Series(list(config['samples'].keys()))
        s.to_csv(output.samples_file, header=False, index=False)


rule gen_bcftools_genome_wide_indi:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        vcf = temp(os.path.join(OUTBIG, "vcf", "{sample}_{chrom}.bcf.gz"))
    threads: 2
    shell:
        """{BCFTOOLS} mpileup -r {wildcards.chrom}  -B -Q {MIN_BQ}  -q {MIN_MQ} --threads {threads} -O u --fasta-ref {REF} --per-sample-mF -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR {input}  | {BCFTOOLS} call -Ob -o {output} --threads {threads} -c
#{BCFTOOLS} index {output}
"""


#rule index_genome_wide_indi:
#    input:
#        vcf = rules.gen_bcftools_genome_wide_indi.output.vcf
#    output:
#        csi = rules.gen_bcftools_genome_wide_indi.output.vcf + ".csi"
#    shell: "{BCFTOOLS} index {input.vcf}"


rule merge_indi_per_chrom:
    input:
        vcf = expand(os.path.join(OUTBIG, "vcf", "{sample}_{{chrom}}.bcf.gz"), sample = config["samples"].keys()),
        csi = expand(os.path.join(OUTBIG, "vcf", "{sample}_{{chrom}}.bcf.gz.csi"), sample = config["samples"].keys()),
        samples_list = os.path.join(OUTBIG, "vcf", "samples.list")
    output:
        vcf = temp(os.path.join(OUTBIG, "vcf", "merged_{chrom}.bcf.gz"))
    threads: 3
    shell:
        "{BCFTOOLS} merge --threads 3 --force-samples -Ou {input.vcf} | {BCFTOOLS} reheader -s {input.samples_list} | {BCFTOOLS} view --threads 3 -Ob -o {output.vcf}"


rule concat_chroms:
    input:
        expand(os.path.join(OUTBIG, "vcf", "merged_{chrom}.bcf.gz"), chrom = config["chroms"])
    output:
        vcf=os.path.join(OUTBIG, "vcf", "merged_snps.bcf.gz"),
    shell: """
    {BCFTOOLS} concat --naive -Ob -o {output.vcf} {input}
    """


rule index_bcf:
    input:
        vcf=os.path.join(OUTBIG, "vcf", "merged_snps.bcf.gz"),
    output:
        csi=os.path.join(OUTBIG, "vcf", "merged_snps.bcf.gz.csi")
    shell: """
    {BCFTOOLS} index {input.vcf}
    """



rule clean_bcf:
    input:
        vcf = os.path.join(OUTBIG, "vcf", "merged_snps.bcf.gz"),
        csi = os.path.join(OUTBIG, "vcf", "merged_snps.bcf.gz.csi")
    output:
        vcf = os.path.join(OUTBIG, "vcf", "merged_snps_filtered.bcf.gz"),
    params:
        bed = config["filters"]["bed"],
    threads: 3
    shell:"""
    {BCFTOOLS} view --threads {threads} -T {params.bed} --include 'STRLEN(REF)=1 & (STRLEN(ALT)=1 | ALT=".")' --max-alleles 2 -Ob -o {output.vcf} {input.vcf}
    """



rule index_clean_bcf:
    input:
        vcf=os.path.join(OUTBIG, "vcf", "merged_snps_filtered.bcf.gz"),
    output:
        csi=os.path.join(OUTBIG, "vcf", "merged_snps_filtered.bcf.gz.csi")
    shell: """
    {BCFTOOLS} index {input.vcf}
    """




# rules for masking low depth genos and hets with low allele support, copied from malthe
rule mask_low_depth:
    """ Set low-depth genotypes as missing and recompute tags """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_{dp}dp.bcf.gz",
    log:
        "{path}_{dp}dp.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} +setGT \
                --output-type u \
                --threads {threads} \
                {input.vcf} \
                -- \
                --include "FMT/DP<{wildcards.dp}" \
                --target-gt q \
                --new-gt . | \
            {BCFTOOLS} +fill-tags \
                --output-type b \
                --threads {threads} \
                -- \
                --tags AN,AC,MAF \
                > {output.vcf} \
        ) 2> {log}
        """


rule mask_het_allele_support:
    """ Set heterozygous genotypes with low support as missing and recompute tags """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_{het}het.bcf.gz",
    log:
        "{path}_{het}het.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} +setGT \
                --output-type u \
                --threads {threads} \
                {input.vcf} \
                -- \
                --include "FMT/GT=='het' & (FMT/AD[*:0]<{wildcards.het} | FMT/AD[*:1]<{wildcards.het})" \
                --target-gt q \
                --new-gt . | \
            {BCFTOOLS} +fill-tags \
                --output-type b \
                --threads {threads} \
                -- \
                --tags AN,AC,MAF \
                > {output.vcf} \
        ) 2> {log}
        """

        

rule filter_missing:
    """ Filter to retain only sites with no missing genotypes """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_nomissing.bcf.gz",
    log:
        "{path}_nomissing.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --exclude 'GT[*] = "mis"' \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
            > {output.vcf} \
        ) 2> {log}
        """


rule filter_variable:
    """ Filter to retain only variable sites. """
    input:
        vcf="{path}.bcf.gz",
    output:
        vcf="{path}_variable.bcf.gz",
    log:
        "{path}_variable.log",
    threads: 10
    shell:
        """
        ( \
            {BCFTOOLS} view \
                --min-alleles 2 \
                --output-type b \
                --threads {threads} \
                {input.vcf} \
                > {output.vcf} \
        ) 2> {log}
        """


rule index:
    """ Index VCF """
    input:
        vcf="{path}.bcf.gz",
    output:
        csi="{path}.bcf.gz.csi",
    threads: 4
    shell:
        """
        {BCFTOOLS} index --threads {threads} {input.vcf}
        """

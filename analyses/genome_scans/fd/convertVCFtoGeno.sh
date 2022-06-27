
BCFTOOLS=/home/genis/software/bcftools/bcftools
CONVERT=/home/genis/software/genomics_general/VCF_processing/parseVCF.py

vcf=/home/genis/warthogs/analyses/dstats_pig/windows_dstats/warts_highdepth_all.vcf.gz
geno=warts_highdepth_all.geno.gz

python3 $CONVERT -i $vcf | bgzip > $geno


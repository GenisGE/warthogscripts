#https://github.com/simonhmartin/genomics_general#abba-baba-statistics-in-sliding-windows
BCFTOOLS=/home/genis/software/bcftools/bcftools
CONVERT=/home/genis/software/genomics_general/VCF_processing/parseVCF.py
DWIN=/home/genis/software/genomics_general/ABBABABAwindows.py


#bcf=/steveData/genis/warthogs/geno_calls_forqpgraph/vcf/merged_snps_filtered.bcf.gz
geno=/home/genis/warthogs/analyses/dstats_pig/windows2/warts_highdepth_all.geno.gz
pops=/home/genis/warthogs/analyses/genome_scan/fd_basedv231052021/pops_combined.list

p2=SE
out=abbababa_win100kb_Ghana${p2}Desert_combined.csv.gz
python3 $DWIN -g $geno --windType coordinate --windSize 100000 --minData 1 -P1 Ghana -P2 $p2 -P3 Desert -O DomesticPig --popsFile $pops --addWindowID --writeFailedWindows -f phased -o $out

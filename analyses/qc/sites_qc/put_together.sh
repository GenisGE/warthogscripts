BEDTOOLS=/home/genis/software/bedtools2/bin/bedtools
ANGSD=/home/genis/software/angsd/angsd


auto=/davidData/genis/references/pig/susscrofa_autosomes.bed
map=/davidData/genis/references/pig/mappability/results/bed/mappability_m1_k100_e2.bed
rep=/home/jfl323/paper/pigMappedNonRepeat.bed
het=/home/jfl323/paper/sites/results/keeplistInbreedSites_minF-0.9_rm20000k.BED
dep=/home/genis/warthogs/genome_masks/depth/results/bed/all_keep.bed


mkdir beds
cd beds

$BEDTOOLS intersect -a $auto -b $rep > autosome_rep.bed
$BEDTOOLS intersect -a autosome_rep.bed -b $het > autosome_rep_het.bed
$BEDTOOLS intersect -a autosome_rep_het.bed -b $dep > autosome_rep_het_dep.bed
$BEDTOOLS intersect -a autosome_rep_het_dep.bed -b $map > autosome_rep_het_dep_map.bed



awk '{print $1"\t"$2+1"\t"$3}'  autosome_rep_het_dep_map.bed > autosome_rep_het_dep_map.regions
$ANGSD sites index  autosome_rep_het_dep_map.regions


# do list that includes X chromosome too
chroms=/davidData/genis/references/pig/susscrofa_chromosomes.bed



$BEDTOOLS intersect -a $chroms -b $rep > chroms_rep.bed
$BEDTOOLS intersect -a chroms_rep.bed -b $map > chroms_rep_map.bed


awk '{print $1"\t"$2+1"\t"$3}' chroms_rep_map.bed > chroms_rep_map.regions
$ANGSD sites index chroms_rep_map.regions



# get mask of removed sites

$BEDTOOLS subtract -a $auto -b autosome_rep_het_dep_map.bed > mask.bed

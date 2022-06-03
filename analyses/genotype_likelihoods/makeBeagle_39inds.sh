ANGSD=/home/genis/software/angsd/angsd

dir=/home/jfl323/paper/beagle/forAnalyses
chrs=/home/jfl323/paper/refAutosomes.rf
bgl=$dir/39samplesNoRepDepthMapHWE_filtered
bams=/home/jfl323/paper/filelists/allFromThesis.filelist
sites=/home/genis/warthogs/genome_masks/all/beds/autosome_rep_het_dep_map.regions

nscaff=`wc -l $chrs | cut -f1 -d" "`


for n in `seq 1 $nscaff`
do
    chr=`head -$n $chrs | tail -1`
    bgl=tmpbgl_$n
    echo "$ANGSD -GL 2 -out $bgl -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -sites $sites -bam $bams -minmapQ 30 -minQ 30 -r $chr -minmaf 0.05"
done | xargs -L1 -I CMD -P20 nice bash -c CMD


echo "Finished estiamting genotype likelihoods going to concatenate all scaffolds in $bgl.beagle"

zcat tmpbgl_1.beagle.gz | head -1 > $bgl.beagle
for i in `seq 1 $nscaff`
do
    zcat tmpbgl_$i.beagle.gz | sed '1d' >> $bgl.beagle
    cat tmpbgl_${i}.arg >> $bgl.arg
    echo "done scaffold $i"
done

gzip $bgl.beagle

rm tmpbgl_*

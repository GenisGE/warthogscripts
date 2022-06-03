EVALADMIX=/home/genis/software/evalAdmix/evalAdmix

bgl=/home/jfl323/paper/beagle/forAnalyses/35samplesNoRepDepthMapHWE_filtered.beagle.gz

resdir=/home/genis/warthogs/analyses/ngsadmix/results
outdir=/home/genis/warthogs/analyses/ngsadmix/evaladmix
outname=wartsOnlycommon35evalAdmix

for k in `seq 2 5`
do
    qfile=`find $resdir/$k | grep qopt_conv`
    ffile=`find $resdir/$k | grep fopt_conv`
    out=$outdir/${outname}_K${k}
    $EVALADMIX -beagle $bgl -qname $qfile -fname $ffile -P 20 -o $out.corres &> $out.log
done

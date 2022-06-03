input=/home/users/xi/Pig_project/fastsimcoal2/1.49e_8_per_generation

for i in {1..50}
do

mkdir $input/$i
cp $input/admixture1.est $input/$i/$i.est
cp $input/admixture1.tpl $input/$i/$i.tpl
cp $input/admixture1_jointMAFpop1_0.obs $input/$i/${i}_jointMAFpop1_0.obs
cp $input/admixture1_jointMAFpop2_0.obs $input/$i/${i}_jointMAFpop2_0.obs
cp $input/admixture1_jointMAFpop3_0.obs $input/$i/${i}_jointMAFpop3_0.obs
cp $input/admixture1_jointMAFpop2_1.obs $input/$i/${i}_jointMAFpop2_1.obs
cp $input/admixture1_jointMAFpop3_1.obs $input/$i/${i}_jointMAFpop3_1.obs
cp $input/admixture1_jointMAFpop3_2.obs $input/$i/${i}_jointMAFpop3_2.obs

export DIR=$input/$i
cd $DIR

/home/users/xi/software/fsc26_linux64/fsc26 -t $i.tpl -n100000 -e $i.est -M -L20 -q -m -c10

rm $i.est
rm $i.tpl
rm ${i}_jointMAFpop1_0.obs
rm ${i}_jointMAFpop2_0.obs
rm ${i}_jointMAFpop3_0.obs
rm ${i}_jointMAFpop2_1.obs
rm ${i}_jointMAFpop3_1.obs
rm ${i}_jointMAFpop3_2.obs

done

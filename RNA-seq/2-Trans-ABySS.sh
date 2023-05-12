#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J transabyssBBD
#SBATCH --mem 32768
#SBATCH -o slurm/transabyssBBD.%j.out
#SBATCH -e slurm/transabyssBBD.%j.err 

kmer1=21
kmer2=29
kmer3=39
kmer4=59

cd /mnt/lustre/home/belben01/transabyss-master/non-rRNA
for read1 in $(ls | grep '_R1'); 
do
	echo $read1
	sample_name=${read1%%_*}; 
	echo $sample_name;
	name1=test.k${kmer1}_${sample_name}
	name2=test.k${kmer2}_${sample_name}
	name3=test.k${kmer3}_${sample_name}
	name4=test.k${kmer4}_${sample_name}
	
	assemblydir1=./${sample_name}/${name1}
	assemblydir2=./${sample_name}/${name2}
	assemblydir3=./${sample_name}/${name3}
	assemblydir4=./${sample_name}/${name4}
	finalassembly1=${assemblydir1}/${name1}-final.fa
	finalassembly2=${assemblydir2}/${name2}-final.fa
	finalassembly3=${assemblydir3}/${name3}-final.fa
	finalassembly4=${assemblydir4}/${name4}-final.fa
	
	mergedassembly=./${sample_name}/${sample_name}_merged.fa

	read2=$(echo $read1| sed 's/_R1/_R2/g');
	echo $read2;


	/mnt/lustre/home/belben01/transabyss-master/transabyss -k ${kmer1} --se ${read1} ${read2} --outdir ${assemblydir1} --name ${name1} --threads 2;

	/mnt/lustre/home/belben01/transabyss-master/transabyss -k ${kmer2} --se ${read1} ${read2} --outdir ${assemblydir2} --name ${name2} --threads 2;
	
	/mnt/lustre/home/belben01/transabyss-master/transabyss -k ${kmer3} --se ${read1} ${read2} --outdir ${assemblydir3} --name ${name3} --threads 2;
	
	/mnt/lustre/home/belben01/transabyss-master/transabyss -k ${kmer4} --se ${read1} ${read2} --outdir ${assemblydir4} --name ${name4} --threads 2;

	/mnt/lustre/home/belben01/transabyss-master/transabyss-merge --mink ${kmer1} --maxk ${kmer4} --prefixes k${kmer1}. k${kmer2}. k${kmer3}. k${kmer4}. --out ${mergedassembly} ${finalassembly1} ${finalassembly2} ${finalassembly3} ${finalassembly4};

done;

#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J assembly-stats
#SBATCH --mem 4096
#SBATCH -o slurm/stats.%j.out
#SBATCH -e slurm/stats.%j.err

# copio los merged fasta del transabyss en otra carpeta para hacer los assembly-stats
cd /mnt/lustre/home/belben01/transabyss-master/non-rRNA
for folder in *
	do
		for file in $folder/*_merged.fa
			do
				cp $file /mnt/lustre/home/belben01/assembly-stats/
			done
	done

cd /mnt/lustre/home/belben01/transabyss_merged_fasta

for file in $(ls | grep '.fa')
	do
		sample_name=${file%%.fa*}
		assembly-stats -t $file > /mnt/lustre/home/belben01/ass-stats/stats_$sample_name.txt
	done

file1=(*)
head -1 $file1 > headers.txt #store the header of the first document (it does not matter which one) from directory	
tail -q -n +2 *.txt > /mnt/lustre/home/belben01/ass-stats/merged_assembly.txt #merge all content from files removing headers (only picks second row)
cat headers.txt merged_assembly.txt > final_assembly.txt #final table with all samples and headers

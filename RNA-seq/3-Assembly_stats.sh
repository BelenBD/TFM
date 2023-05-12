#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J assembly-stats
#SBATCH --mem 4096
#SBATCH -o slurm/stats.%j.out
#SBATCH -e slurm/stats.%j.err

# Copy merged fasta files from Transabyss to another directory for assembly-stats
cd /mnt/lustre/home/belben01/transabyss-master/non-rRNA
for folder in *
	do
		for file in $folder/*_merged.fa
			do
				cp $file /mnt/lustre/home/belben01/assembly-stats/
			done
	done

# Change directory to the copied merged fasta files location
cd /mnt/lustre/home/belben01/transabyss_merged_fasta

# For each file with ".fa" extension in the current directory, run assembly-stats and save the output in a file with sample name in another directory
for file in $(ls | grep '.fa')
	do
		sample_name=${file%%.fa*}
		assembly-stats -t $file > /mnt/lustre/home/belben01/ass-stats/stats_$sample_name.txt
	done

# Store the header of the first document (it does not matter which one) from directory	
file1=(*)
head -1 $file1 > headers.txt

# Merge all content from files removing headers (only picks second row) and save it to a file with merged assembly statistics
tail -q -n +2 *.txt > /mnt/lustre/home/belben01/ass-stats/merged_assembly.txt

# Merge the headers file and merged assembly statistics file to create a final table with all samples and headers
cat headers.txt merged_assembly.txt > final_assembly.txt 

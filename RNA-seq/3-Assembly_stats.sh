#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J assembly-stats
#SBATCH --mem 4096
#SBATCH -o slurm/stats.%j.out
#SBATCH -e slurm/stats.%j.err

# Copy merged fasta files from Transabyss to another directory for assembly-stats
find /mnt/lustre/home/belben01/transabyss-master/non-rRNA -name "*_merged.fa" -exec cp {} /mnt/lustre/home/belben01/assembly-stats/ \;

# Change directory to the copied merged fasta files location
cd /mnt/lustre/home/belben01/assembly-stats

# For each file with ".fa" extension in the current directory, run assembly-stats and save the output in a file with sample name in another directory
for file in *.fa
	do
		sample_name=${file%%.fa*}
		assembly-stats -t $file > /mnt/lustre/home/belben01/ass-stats/stats_$sample_name.txt
	done

# Store the header of the first document (it does not matter which one) from directory	
head -1 /mnt/lustre/home/belben01/ass-stats/stats_*.txt > headers.txt

# Merge all content from files removing headers (only picks second row) and save it to a file with merged assembly statistics
tail -q -n +2 /mnt/lustre/home/belben01/ass-stats/stats_*.txt > /mnt/lustre/home/belben01/ass-stats/merged_assembly.txt

# Merge the headers file and merged assembly statistics file to create a final table with all samples and headers
printf '%s\n' "$(cat headers.txt)" "$(cat merged_assembly.txt)" > final_assembly.txt


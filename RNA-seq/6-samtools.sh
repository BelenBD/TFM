#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J samtools
#SBATCH --mem 52768
#SBATCH -t 7-02:30:15 
#SBATCH -o slurm/samtools.%j.out
#SBATCH -e slurm/samtools.%j.err 


cd /mnt/lustre/home/belben01/bwa-sam

for file in $(ls)
	do
		input=${file}
		echo $file
		sample_name=${file%%.sam*}
		echo $sample_name
		samtools view -b -S /mnt/lustre/home/belben01/bwa-sam/$file > /mnt/lustre/home/belben01/bwa-bam/$sample_name.bam
    		#samtools sort /mnt/lustre/home/belben01/bwa-bam/$sample_name.bam -o /mnt/lustre/home/belben01/bwa-bam/Sorted/$sample_name-sorted.bam
	done

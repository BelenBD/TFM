#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J count
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/count.%j.out
#SBATCH -e slurm/count.%j.err

cd /mnt/lustre/home/belben01/bwa-bam

for file in $(ls| grep '.bam')
	do
		echo $file
		sample_name=${file%%.bam*}
		gffread /mnt/lustre/home/belben01/eggnogmapper/$sample_name/$sample_name.emapper.gff -F -T -o /mnt/lustre/home/belben01/bwa-bam/gffread/$sample_name.gtf
	done


cd /mnt/lustre/home/belben01/bwa-bam

for file in $(ls | grep '.bam')
	do
		bam_file=${file}
		echo $bam_file
		sample_name=${file%%.bam*}
		echo $sample_name
		featureCounts -a /mnt/lustre/home/belben01/bwa-bam/gffread/$sample_name.gtf -o /mnt/lustre/home/belben01/featurecounts/$sample_name.txt $bam_file -t transcript -g em_target
	done

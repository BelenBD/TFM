#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J bwa
#SBATCH --mem 32768
#SBATCH -o slurm/bwa.%j.out
#SBATCH -e slurm/bwa.%j.err 

cd /mnt/lustre/home/belben01/eggnogmapper

for folder in $(ls)
	do
		sample_name=${folder}
		R1=${sample_name}_R1.fq.gz
		R2=${sample_name}_R2.fq.gz
		mkdir /mnt/lustre/home/belben01/index-bwa/$sample_name
		bwa index -p /mnt/lustre/home/belben01/index-bwa/$sample_name/$sample_name -a is /mnt/lustre/home/belben01/STAR/Samples/$sample_name.fa
		bwa mem /mnt/lustre/home/belben01/index-bwa/$sample_name/$sample_name /mnt/lustre/home/belben01/sortmerna/Reads/non-rRNA/$R1 /mnt/lustre/home/belben01/sortmerna/Reads/non-rRNA/$R2 > /mnt/lustre/home/belben01/bwa-sam/$sample_name.sam

	done
	

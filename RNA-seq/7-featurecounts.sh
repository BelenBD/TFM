#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J count
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/count.%j.out
#SBATCH -e slurm/count.%j.err

cd /mnt/lustre/home/belben01/bwa-bam

find . -name '*.bam' | while read -r bam_file; do
    echo "$bam_file"
    sample_name="${bam_file%%.bam*}"
    echo "$sample_name"
    gffread /mnt/lustre/home/belben01/eggnogmapper/"$sample_name"/"$sample_name".emapper.gff -F -T -o gffread/"$sample_name".gtf
    featureCounts -a gffread/"$sample_name".gtf -o /mnt/lustre/home/belben01/featurecounts/"$sample_name".txt "$bam_file" -t transcript -g em_target
done

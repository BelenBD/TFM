#!/bin/bash
#SBATCH -p short
#SBATCH -n 16
#SBATCH -J fastqc
#SBATCH --mem 32768  
#SBATCH -o slurm/fastqc.%j.out
#SBATCH -e slurm/fastqc.%j.err

INPUT_DIR="/mnt/lustre/home/belben01/Raw"
OUTPUT_DIR="/mnt/lustre/home/belben01/FastQC"

# Loop over all sequenced reads from /Raw folder

for file in $INPUT_DIR/*R1*; do
  file_R2=${file/R1/R2}

  # Executing FastQC in both reads for each sample
  fastqc $file $file_R2 -o $OUTPUT_DIR
done

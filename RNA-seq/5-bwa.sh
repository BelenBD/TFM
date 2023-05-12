#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J bwa
#SBATCH --mem 32768
#SBATCH -o slurm/bwa.%j.out
#SBATCH -e slurm/bwa.%j.err 

# Define input and output directories
input_dir="/mnt/lustre/home/belben01/eggnogmapper"
index_dir="/mnt/lustre/home/belben01/index-bwa"
sam_dir="/mnt/lustre/home/belben01/bwa-sam"

# Iterate over input directories
for folder in "${input_dir}"/*
do
    if [[ -d "${folder}" ]]; then
        sample_name="$(basename "${folder}")"
        R1="${sample_name}_R1.fq.gz"
        R2="${sample_name}_R2.fq.gz"
        index_file="${index_dir}/${sample_name}/${sample_name}"
        sam_file="${sam_dir}/${sample_name}.sam"

        # Create index directory
        mkdir -p "${index_dir}/${sample_name}"

        # Index reference file with BWA
        printf "Indexing reference file for %s...\n" "${sample_name}"
        bwa index -p "${index_file}" -a is "${input_dir}/${sample_name}.fa"

        # Map reads with BWA
        printf "Mapping reads for %s...\n" "${sample_name}"
        bwa mem "${index_file}" \
            "${input_dir}/${R1}" \
            "${input_dir}/${R2}" \
            > "${sam_file}"
    fi
done

echo "Done mapping reads with BWA."

	

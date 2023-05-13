#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J samtools
#SBATCH --mem 52768
#SBATCH -t 7-02:30:15 
#SBATCH -o slurm/samtools.%j.out
#SBATCH -e slurm/samtools.%j.err 

input_dir="/mnt/lustre/home/belben01/bwa-sam"
output_dir="/mnt/lustre/home/belben01/bwa-bam"

# Check if output directories exist, if not create them
for dir in "$output_dir" "$output_dir/Sorted"; do
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
  fi
done

# Loop through sam files in input directory
for file in "$input_dir"/*.sam; do
  [[ -f "$file" ]] || continue  # skip if not a file
  echo "Processing $file"

  # Extract sample name from file name
  sample_name="$(basename "${file%.sam}")"
  echo "Sample name: $sample_name"

  # Convert sam to bam
  samtools view -b "$file" > "$output_dir/$sample_name.bam"

  # Sort bam file
  samtools sort "$output_dir/$sample_name.bam" -o "$output_dir/Sorted/$sample_name-sorted.bam"
done


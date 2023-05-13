#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J count
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/count.%j.out
#SBATCH -e slurm/count.%j.err

# Define variables for directories
BAM_DIR="/mnt/lustre/home/belben01/bwa-bam"
GFFREAD_DIR="/mnt/lustre/home/belben01/gffread"
FEATURECOUNTS_DIR="/mnt/lustre/home/belben01/featurecounts"
EGGNOGMAPPER_DIR="/mnt/lustre/home/belben01/eggnogmapper"

# Change to the BAM directory
cd "$BAM_DIR"

# Loop through all BAM files in the directory and its subdirectories
find . -name '*.bam' -type f | while read -r bam_file; do
    # Extract sample name from the file name
    sample_name="${bam_file%.bam}"
    sample_name="${sample_name##*/}"
    echo "Processing sample: $sample_name"
    
    # Run gffread to convert the eMapper output GFF file to GTF format
    gffread "$EGGNOGMAPPER_DIR/$sample_name/$sample_name.emapper.gff" \
        -F -T -o "$GFFREAD_DIR/$sample_name.gtf"
    
    # Run featureCounts to count reads that overlap with transcripts
    featureCounts -a "$GFFREAD_DIR/$sample_name.gtf" \
        -o "$FEATURECOUNTS_DIR/$sample_name.txt" \
        "$bam_file" -t transcript -g em_target
done


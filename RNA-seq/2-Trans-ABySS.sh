#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J transabyssBBD
#SBATCH --mem 32768
#SBATCH -o slurm/transabyssBBD.%j.out
#SBATCH -e slurm/transabyssBBD.%j.err 

#!/bin/bash

# Set kmer sizes
kmers=(21 29 39 59)

# Set input and output directories
input_dir=/mnt/lustre/home/belben01/transabyss-master/non-rRNA
output_dir=/mnt/lustre/home/belben01/transabyss-master/non-rRNA/merged

# Create output directory if it does not exist
mkdir -p $output_dir

# Loop through all read1 files in input directory
for read1 in $input_dir/*_R1.fastq.gz; do
    echo "Processing $read1..."

    # Extract sample name from read1 file name
    sample_name=$(basename $read1 | cut -d '_' -f 1)

    # Loop through all kmer sizes
    for kmer in ${kmers[@]}; do
        echo "Running transabyss with kmer size $kmer..."

        # Set output directory and file names
        assemblydir=$output_dir/${sample_name}_k${kmer}
        finalassembly=${assemblydir}/final.fa

        # Run transabyss with current kmer size
        transabyss --se $read1 ${read1/_R1/_R2} --outdir $assemblydir --name ${sample_name}_k${kmer} --threads 2 -k $kmer

        # Check if transabyss produced output file
        if [ ! -f $finalassembly ]; then
            echo "Error: transabyss did not produce output file"
            exit 1
        fi
    done

    # Set merged assembly file name
    mergedassembly=$output_dir/${sample_name}_merged.fa

    # Merge assemblies using transabyss-merge
    transabyss-merge --mink ${kmers[0]} --maxk ${kmers[-1]} --prefixes ${kmers[@]/#/k} --out $mergedassembly ${output_dir}/${sample_name}_k*/final.fa
done


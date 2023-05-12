#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J buscoBBD
#SBATCH --mem 32768
#SBATCH -o slurm/buscoBBD.%j.out
#SBATCH -e slurm/buscoBBD.%j.err 

# Copy the merged fasta files from transabyss to another directory for BUSCO analysis
cd /mnt/lustre/home/belben01/transabyss-master/non-rRNA
for folder in * 
  do
    for file in $folder/*_merged.fa
      do
        cp $file /mnt/lustre/home/belben01/BUSCO/
    done
done

# Go to the directory where the fasta files were copied
cd /mnt/lustre/home/belben01/transabyss_merged_fasta

# Run BUSCO on each fasta file
for file in $(ls | grep '.fa')
	do
		input=${file}
		output=${file%%.*} 
		busco -i $input -l fungi_odb10 -m tran -o $output -c 12 # run BUSCO on the input file with options for the fungi database, transcript mode, and number of threads
	done


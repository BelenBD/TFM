#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J buscoBBD
#SBATCH --mem 32768
#SBATCH -o slurm/buscoBBD.%j.out
#SBATCH -e slurm/buscoBBD.%j.err 

# copio los merged fasta del transabyss en otra carpeta para hacer el análisis BUSCO
cd /mnt/lustre/home/belben01/transabyss-master/non-rRNA
for folder in *
  do
    for file in $folder/*_merged.fa
      do
        cp $file /mnt/lustre/home/belben01/BUSCO/
    done
done

cd /mnt/lustre/home/belben01/transabyss_merged_fasta

for file in $(ls | grep '.fa')
	do
		input=${file}
		output=${file%%.*}
		busco -i $input -l fungi_odb10 -m tran -o $output -c 12	
	done
	

#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J emapper
#SBATCH --mem 12288
#SBATCH -o slurm/emapper.%j.out
#SBATCH -e slurm/emapper.%j.err 

# Download the necessary data for EggNOG mapper
download_eggnog_data.py

# Create a diamond database for the Fungi taxonomic group
create_dbs.py -m diamond --dbname fungi --taxa Fungi

# For each input FASTA file in the transabyss_merged_fasta directory, run EggNOG mapper with the diamond database and other options
cd /mnt/lustre/home/belben01/transabyss_merged_fasta
for file in $(ls | grep '.fa')
	do
		input=${file} 
		sample_name=${file%%.fa*} 
		folder=${sample_name} 
		mkdir /mnt/lustre/home/belben01/eggnogmapper/$folder 
		emapper.py -m diamond --itype CDS -i $input -o /mnt/lustre/home/belben01/eggnogmapper/$folder/$sample_name --dmnd_db /mnt/lustre/home/belben01/anaconda3/envs/eggnogmapper/lib/python3.7/site-packages/data/fungi.dmnd --tax_scope Saccharomycetaceae --decorate_gff yes --cpu 0
	done

	

#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J rRNA-LM
#SBATCH --mem 32768  
#SBATCH -o slurm/sortmerna_LM.%j.out
#SBATCH -e slurm/sortmerna_LM.%j.err

#Indexed databases:
indexdb_rna --ref rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:\
rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:\
rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db

for read1 in Raw/*_R1*; 
do read2=$(echo $read1| sed s/'_R1_'/'_R2_'/);  
sample_name=$(echo $read1| grep -oP '(?<=/).*(?=_S)'); 
echo $sample_name;  
echo $read1;  
echo $read2;

sortmerna --ref Metatranscriptomica/silva-euk-18s-id95.fasta --ref Metatranscriptomica/silva-euk-28s-id98.fasta --ref Metatranscriptomica/rfam-5.8s-database-id98.fasta \
--reads $read1 --reads $read2 --out2 --other --fastx --paired_in \
--workdir sortmerna/;


mv sortmerna/out sortmerna/$sample_name;

mv sortmerna/$sample_name/other_fwd.fq.gz sortmerna/Reads/non-rRNA/$sample_name"_R1.fq.gz";
mv sortmerna/$sample_name/other_rev.fq.gz sortmerna/Reads/non-rRNA/$sample_name"_R2.fq.gz";

rm -rf sortmerna/kvdb;

done

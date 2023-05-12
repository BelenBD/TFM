#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J kraken
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/kraken.%j.out
#SBATCH -e slurm/kraken.%j.err

# Run Kraken2 on paired-end reads
for read1 in reads/*R1_001_72nt.fastq.gz; do 
    read2=$(echo $read1| sed s/'_R1_'/'_R2_'/); 
    sample_name=$(echo $read1| grep -oP '(?<=/).*(?=_S)'); 
    echo "Running Kraken2 on $sample_name"; 
    $HOME/kraken2/kraken2 --use-names --threads 16 --db Kraken2/Wine/vino --report kraken2_out/$sample_name.report.txt --paired $read1 $read2 > kraken2_out/$sample_name.kraken;
done

# Run Bracken on Kraken2 report files
for report in kraken2_out/*.report.txt; do 
    outp=$(echo $report| grep -oP '(?<=/).*(?=.report)'); 
    echo "Running Bracken on $outp"; 
    $HOME/Bracken-2.5/bracken -d Kraken2/Wine/vino -i $report -r 72 -l S -o bracken_out/$outp.bracken;
done

# Filter and concatenate Bracken output files
mkdir -p filtered_bracken_out
for file in bracken_out/*.bracken; do
    name=$(basename "$file" .bracken)
    awk -v name="$name" '$4 == "S" && $1 > 0 {print $0, name}' "$file" > "filtered_bracken_out/${name}_filtered.txt"
done

cd filtered_bracken_out
awk '{print $1,$2,$3,$4,$5,$6,FILENAME}' *_filtered.txt | sed 's/_filtered.txt//' > all_filtered_bracken.txt

#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=domout
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1080
#SBATCH --time=07:00:00
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger

module purge
module add HMMER/3.3.2-gompi-2020b

input_list="../Genome_to_species.txt"

readarray list <$input_list
file=${list[${SLURM_ARRAY_TASK_ID} - 1]}
spe=$(echo $file | cut -d "/" -f1)
gen=$(echo $file | cut -d "/" -f2)

mkdir output
mkdir output/domtables
mkdir output/domtables/domtable_${spe}
for fam in $(cat input/lists/gene_subfamily_list.txt); do 
	echo ${fam}_v
	hmmsearch --domtblout output/domtables/domtable_${spe}/${fam}.${spe}.domout.tab --noali --notextw -E 1e-3 --domE 1e-3 --incE 0.001 input/gene_subfamilies/${fam}/${fam}.ginsi.hmm input/protein/${spe}/${spe}_protein.faa
done
echo "domout done"

#!/bin/bash

#SBATCH --job-name=6ft
#SBATCH --mail-type=none
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=350G
#SBATCH --time=03:30:00
#SBATCH --qos=standard
#SBATCH --partition=scavenger,main

module purge
module add Biopython/1.79-foss-2021a

input_list="../Genome_to_species.txt"

readarray list <$input_list
file=${list[${SLURM_ARRAY_TASK_ID} - 1]}
spe=$(echo $file | cut -d "/" -f1)
gen=$(echo $file | cut -d "/" -f2)


mkdir ./input/protein/${spe}
python3 input/python_codes/sft_beta.py -f ${spe}_genome.fna -p 15 -o input/protein/${spe}/ -i input/target_genomes

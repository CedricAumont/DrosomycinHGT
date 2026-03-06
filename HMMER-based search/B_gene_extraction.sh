#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=Job2-3-extraction.sh
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G
#SBATCH --time=00:30:00
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger

module purge
module add Biopython/1.79-foss-2022a

input_list="input/lists/species_list.txt"

readarray list <$input_list
SPE2=${list[${SLURM_ARRAY_TASK_ID} - 1]}
SPE=$(basename $SPE2)

python3 input/python_codes/gene_extraction_trinity.py -i output/targetcleaned/${SPE} -s input/protein/${SPE} -g /scratch/ca3321fu/drosomycin_project/TSA_assemblies/assemblies_blattodea_he/${SPE}_transc.fna -o input/gene_families_2
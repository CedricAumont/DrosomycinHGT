#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=Job2-3-extraction.sh
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50
#SBATCH --time=00:05:00
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger

module purge
module add Biopython/1.79-foss-2022a

input_list="input/lists/species_list.txt"

readarray list <$input_list
SPE2=${list[${SLURM_ARRAY_TASK_ID} - 1]}
SPE=$(basename $SPE2)

mkdir output/targetcleaned_2
python3 input/python_codes/clean_domtable.py -d output/domtables_2/domtable_${SPE} -c 50 -o output/targetcleaned_2
#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=extraction
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=510G
#SBATCH --time=02:00:00
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger

module purge
module add Biopython/1.79-foss-2022a
input_list="../species_w_drs2.txt"

readarray list <$input_list
file=${list[${SLURM_ARRAY_TASK_ID} - 1]}
spe=$(basename $file)

echo ${spe}

python3 input/python_codes/gene_extraction_only.py -i output/targetcleaned_2/${spe} -s input/protein/${spe} -g input/target_genomes/${spe}_genome.fna -o input/gene_subfamilies_out

echo "done"
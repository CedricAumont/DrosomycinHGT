#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=J7
#SBATCH --mail-type=none
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --time=00:30:00
#SBATCH --qos=standard
#SBATCH --partition=main,scavenger

module purge
module add HMMER/3.3.2-gompi-2020b
module add MAFFT

input_list="input/lists/gene_subfamily_list.txt"

readarray list <$input_list
FAM2=${list[${SLURM_ARRAY_TASK_ID} - 1]}
FAM=$(basename $FAM2)

mkdir input/gene_subfamilies_out
mkdir input/gene_subfamilies_out/${FAM}

mafft --globalpair --maxiterate 1000 input/gene_subfamilies_2/${FAM}/${FAM}.fa >input/gene_subfamilies_2/${FAM}/${FAM}.ginsi.alin
esl-reformat -o input/gene_subfamilies_2/${FAM}/${FAM}.ginsi.sto stockholm input/gene_subfamilies_2/${FAM}/${FAM}.ginsi.alin
echo "alignment done"

hmmbuild input/gene_subfamilies_2/${FAM}/${FAM}.ginsi.hmm input/gene_subfamilies_2/${FAM}/${FAM}.ginsi.sto
echo "hmm profile done"


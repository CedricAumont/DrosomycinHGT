Explanation of the pipeline for genomes.
All jobs are run with Slurm Workload Manager

1) Six frame translation

All genomes are gathered in a single folder (target_genomes) and file name were formatted as such {genus}_{species}_genome.fna
The script Six_frame_translation.sh was performed as a job-array on all genomes. It uses the python script sft_beta.py. 
Proteome outputs are directed in species specific folder 
protein/${spe}/

2) HMM profile building and search

Script J3-msa_and_profile.sh is used to make the HMM profile using the 
list of immune genes (gene_subfamily_list.txt) and the 
orthologuous genes in the immune gene sub-family folder (i.e. Drosomycin genes). The msa, profile and stokholm format alignment is 
provided in the gene-subfamily folder.

Script J4-domout_making.sh perform the HMM-search on every genomes.

Script J5-cleaning.sh parse the output of the HMM-search and use the
python script clean_domtable.py to extract tabular information on the genes making at least 75% of the profile length. It results in a "cleaned" HMM output.

3) Gene extraction

Script J6-extraction.sh launches the python script gene_extraction_only.py to extract the genes selected in the "cleaned" HMM output for all species where a Drosomycin was found. the protein sequences are retrieved in the folder gene_subfamilies_2.

4) Second iteration

To improve the search we perform the step 2) and 3) again using J7, J8, J9, J10. The J9-cleaning2.sh extract tabular information on the genes making at least 50% of the profile length this time. Manual curation of the gene found was performed to remove genes that were not Drs.


Explanation of the pipeline for transcriptome.

The procedure was similar with 2 minor changes:
In step 1) transcriptomes file names were formatted as such {genus}_{species}_genome.fna
in step 3) and 4) the scripts B_gene_extraction.sh and gene_extraction_trinity.py were used.

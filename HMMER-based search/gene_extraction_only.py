# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 12:31:23 2024

@author: caumont
"""


import pandas as pd
import argparse
import os
from collections import defaultdict
import re
from Bio import SeqIO
import ast


def upload_targetclean(fam_name, sp_name, inp):
    tgclean1 = pd.read_csv(f'{inp}/{fam_name}.{sp_name}.targetcleaned.csv', sep = ";")
    with open(f'{inp}/{fam_name}.{sp_name}.geneid.txt', 'r') as fp:
        tgclean2 = ast.literal_eval(fp.read())
    print(tgclean1)
    print(tgclean2)
    return tgclean1, tgclean2


def _parse_sp_fam_name(ofile_path: str) -> tuple[str, str]:
    """
    Parse the species and gene family name for a given file.
    """
    #first_name = ofile_path.strip().split("/")[-1]
    first_name = re.split(r"/|\\", ofile_path.strip())[-1]
    names = first_name.strip().split(".")
    # sp_name, fam_name
    return names[1], names[0]


def get_unique_fa_id(record) -> str:
    """
    creates a unique fasta id for the SeqIO.to_dict function
    """
    oid = record.id # id of the sequence
    strand = re.search(r"strand[FR]", record.description).group(0)
    return f"{oid}_{strand}"


def get_protein_dict(protpath: str) -> dict:
    """
    Create a dictionnary from a fasta file.
    """
    return {
        k: v
        for d in (
            SeqIO.to_dict(
                SeqIO.parse(os.path.join(protpath, fafile), "fasta"),
                key_function=get_unique_fa_id,
            )
            for fafile in os.listdir(protpath)
            if fafile.endswith(".faa")
        )
        for k, v in d.items()
    }


def get_nucl_dict(genome: str) -> dict:
    """
    Create a dictionnary from a fasta file.
    """
    return SeqIO.to_dict(
                SeqIO.parse(genome, "fasta"))


def get_gene(dfandidx, proteomes_dict: dict, file: str, nucl_dict: dict) -> defaultdict[str]:
    sp_name, fam_name = _parse_sp_fam_name(file)
    target_df = dfandidx[0]
    gene_idx = dfandidx[1]
    gene_number = len(gene_idx) // 2
    target_df.index = target_df["index"] # need to remove or replace by level_0 or KEEP
    gene_rec = defaultdict(str)
    nucl_rec = defaultdict(str)
    k = 0

    while k < 2 * gene_number:
        cur_gene_num = k // 2 + 1
        print("Gene ", cur_gene_num)
        idx_st = gene_idx[k]
        idx_end = gene_idx[k + 1]
        strand = target_df.at[idx_st, "strand"]
        scaffold = target_df.at[idx_st, "scaffold"]
        start = target_df.at[idx_st, "pos_st"]
        end = target_df.at[idx_end, "pos_end"]
        descr = target_df.at[idx_st, "description_of_target"].strip().split(", ")[4]
        nucl_header = descr.split(" ")[0]
        print(nucl_header)

        genename = (
            f">{fam_name}_{sp_name}_{cur_gene_num}_{scaffold}_{strand}_{start}_{end}"
        )
        nuclname = (
            f">{fam_name}_{sp_name}_{cur_gene_num}_{scaffold}_{strand}_{start}_{end}"
        )

        for ex in range(idx_st, idx_end + 1):
            exon_name = target_df.at[ex, "target_name"]
            fa_id_match = f"{exon_name}_{strand}" # cf get_unique_fa_id. it looks like >gac123_genus_species_12:154,_strandF
            ali_st = target_df.at[ex, "ALI_from"] - 1 # start of the hit on the exon in aa
            ali_end = target_df.at[ex, "ALI_to"]
            ex_start = target_df.at[ex, "pos_st"] # genomic position of the exon in nt
            ex_end = target_df.at[ex, "pos_end"]
            tar_len = target_df.at[ex, "tlen"]
            trimseq = proteomes_dict[fa_id_match].seq[ali_st:ali_end]
            #print(nucl_dict)
            if strand == "strandR":
                nuclseq = nucl_dict[str(nucl_header)].seq[ex_start+3*((tar_len-ali_end)+1):ex_end-3*ali_st] #scaffolf may not work, depends on the headers of genome
                gene_rec[genename] = trimseq + gene_rec[genename]
                nucl_rec[nuclname] += nuclseq
            else:
                nuclseq = nucl_dict[str(nucl_header)].seq[ex_start+3*ali_st:ex_start+3*ali_end]
                gene_rec[genename] += trimseq
                nucl_rec[nuclname] += nuclseq
        k += 2
    return gene_rec, nucl_rec

def second_fasta(gene_record: defaultdict, outputfile: str):
    with open(outputfile, "a") as opf:
        for gene in gene_record:
            opf.write(f"{gene}\n{gene_record[gene]}\n")


def extraction(infol, pro_path, outfol, nucl_path):
    
    proteomes_dict = get_protein_dict(pro_path)
    nucleo_dict = get_nucl_dict(nucl_path)
    

    file_list = [
        os.path.join(dp, f) for dp, dn, filenames in os.walk(infol) for f in filenames if f.endswith("targetcleaned.csv")
    ]

    for file in file_list:
        sp_name, fam_name = _parse_sp_fam_name(file)
        targetcleaned = upload_targetclean(fam_name, sp_name, infol)

        if len(targetcleaned[0]) > 0:
                    gene_extracted = get_gene(targetcleaned, proteomes_dict, file, nucleo_dict)
                    sec_fafile = f"{outfol}/{fam_name}/{fam_name}.fa"
                    second_fasta(gene_extracted[0], sec_fafile)
                    third_fafile = f"{outfol}/{fam_name}/{fam_name}.fna"
                    second_fasta(gene_extracted[1], third_fafile)
                    
    return 1


def main():
    """
    The main function of the program.
    """
    parser = argparse.ArgumentParser(
        description="Summarises counts per hmm model for all species in a given folder and reconstruct the gene if "
        "asked."
    )
    parser.add_argument(
        "-i",
        "--infol_spe",
        help="The folder in which are the domtout targetcleaned files of the species.",
        type=str,
        default="output/targetcleaned",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--spedict",
        help="the path (only) to the proteome species",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-g",
        "--genome",
        help="the path (incl. filename) to the genome species",
        type=str,
        required=True,
    )
    
    parser.add_argument(
        "-o",
        "--outfol",
        help="The folder name to tsave the amino acid seqences.",
        type=str,
        default="input/gene_families_2",
        required=True,
    )

    
    args = parser.parse_args()

    extraction(
        args.infol_spe,
        args.spedict,
        args.outfol,
        args.genome
    )


if __name__ == "__main__":
    main()
    

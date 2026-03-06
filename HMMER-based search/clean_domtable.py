# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:14:43 2024

@author: caumont
"""

import pandas as pd
import numpy as np
import argparse
import os
from collections import defaultdict
import pickle


def _file_not_empty(annot_file: str) -> bool:
    """
    Return True if the file is not empty.
    """
    if os.path.getsize(annot_file) == 0:
        return False
    with open(annot_file, 'rt') as af:
        for line in af:
            if line.strip() and not line.startswith('#'):
                return True
    return False


def _preproc_description(descr_str: str) -> list:
    """
    Convert a line of hit description into list.
    """
    descr_tmp = descr_str.strip().split(", ")
    # strand + scaffold + start&end pos + description_string
    return (
        [descr_tmp[1]]
        + [descr_tmp[3].split()[1]]
        + descr_tmp[2].split(":")
        + [descr_str]
    )


def sort_domtout(domtout: str) -> pd.DataFrame:
    """
    Sorts the hits of the output of the hmm search per e-value
    and then domain e-value, and converts the file into a pandas dataframe
    """
    match_dict = defaultdict(list)
    column_names = [
        "target_name",
        "taccession",
        "tlen",
        "query_name",
        "qaccession",
        "qlen",
        "FS_Evalue",
        "FS_score",
        "FS_bias",
        "D_nb",
        "D_of",
        "D_c-Evalue",
        "D_i-Evalue",
        "D_score",
        "D_bias",
        "HMM_from",
        "HMM_to",
        "ALI_from",
        "ALI_to",
        "ENV_from",
        "ENV_to",
        "acc",
        "strand",
        "scaffold",
        "pos_st",
        "pos_end",
        "description_of_target",
    ]
    convert_dict = {
        "qlen": int,
        "tlen": int,
        "FS_Evalue": float,
        "D_c-Evalue": float,
        "D_i-Evalue": float,
        "HMM_from": int,
        "HMM_to": int,
        "ENV_from": int,
        "ENV_to": int,
        "ALI_from": int,
        "ALI_to": int,
        "scaffold": int,
        "pos_st": int,
        "pos_end": int,
    }

    with open(domtout, "rt") as domtab:
        for line in domtab:
            if line.strip() and not line.startswith("#"):
                val_row = (spl_line := line.strip().split())[
                    :22
                ] + _preproc_description(" ".join(spl_line[22:]))

                for key, value in zip(column_names, val_row):
                    match_dict[key].append(value)

        result_df = pd.DataFrame(data=match_dict).astype(convert_dict)
        srt_df = result_df.sort_values(
            by=["FS_Evalue", "D_c-Evalue"], ascending=[True, True]
        )
    srt_df.reset_index(inplace=True, drop=True)
    return srt_df


def sort_hits(sorted_df: pd.DataFrame) -> pd.DataFrame:
    """
    Keeps only the best fitting hit to the profile HMM
    and order the hits per scaffold and genomic position.
    """
    L=[]
    for i in range(1,len(sorted_df)):
        if sorted_df["target_name"][i] == sorted_df["target_name"][i-1]:
            L.append(i)
    for i in L:
        sorted_df = sorted_df.drop([i])
    sorted_hits = sorted_df.sort_values(
        by=["scaffold", "pos_st"], ascending=[True, True]
    )
    sorted_hits.reset_index(inplace=True)
    return sorted_hits


def _attribute(targ: int, th_df: pd.DataFrame) -> tuple[int, str, int, int]:
    """
    For a given hit, return scaffold ID, strand, and position on the profile HMM.
    """
    scaff = th_df["scaffold"][targ]
    strand = th_df["strand"][targ]
    pos_st = th_df["HMM_from"][targ]
    pos_end = th_df["HMM_to"][targ]
    return scaff, strand, pos_st, pos_end


def _coverage_test(coverage, percent, qlen, targ, ex_targ, to_keep, to_drop, gene_idx):
    """
    Test if a given sequence covers more than a given percent of the profile HMM.
    """
    if coverage > percent * qlen:
        for i in range(targ, ex_targ + 1):
            to_keep.append(i)
        gene_idx.append(targ)
        gene_idx.append(ex_targ)
    else:
        for i in range(targ, ex_targ + 1):
            to_drop.append(i)


def graph_cleaner(th_df: pd.DataFrame, percent: int) -> tuple[pd.DataFrame, list]:
    """
    graph_cleaner reconstruct the gene-sized hits from the exon-sized hits,
    and remove short hits based on the coverage threshold given. It also keeps
    the index of start and end of each gene-sized hits for gene reconstruction
    """
    to_keep, gene_idx, to_drop = [], [], []
    qlen = th_df["qlen"][0]
    th_df["index"] = [i for i in range(len(th_df))] # This seems to be very important for the pipeline to run properly.

    # case of one target hit
    if len(th_df) == 1:
        coverage_ini = th_df["HMM_to"][0] - th_df["HMM_from"][0]
        _coverage_test(coverage_ini, percent, qlen, 0, 0, to_keep, to_drop, gene_idx)

    # case of several targets
    targ = 0
    while targ < (len(th_df) - 1):
        current_scaffold, current_strand, current_pos_st, current_pos_end = _attribute(
            targ, th_df
        )
        next_scaffold, next_strand, next_pos_st, next_pos_end = _attribute(
            targ + 1, th_df
        )
        coverage = current_pos_end - current_pos_st

        # forward strand
        if current_strand == "strandF":
            if current_scaffold == next_scaffold and current_strand == next_strand:
                # exon loop
                ex_targ = targ
                while (
                    next_pos_st > current_pos_end - (qlen // 2)
                    and current_scaffold == next_scaffold
                    and current_strand == next_strand
                    and ex_targ < (len(th_df) - 1)):
                    ex_targ += 1
                    
                    (current_scaffold,
                        current_strand,
                        current_pos_st,
                        current_pos_end) = _attribute(ex_targ, th_df)
                    if ex_targ != (len(th_df) - 1):
                        (next_scaffold,
                          next_strand, 
                          next_pos_st, 
                          next_pos_end) = _attribute(ex_targ + 1, th_df)
                        coverage += current_pos_end - current_pos_st
                    if ex_targ == (len(th_df) - 1):
                        coverage += current_pos_end - current_pos_st
                _coverage_test(
                    coverage, percent, qlen, targ, ex_targ, to_keep, to_drop, gene_idx
                )
                targ = ex_targ
            else:
                _coverage_test(
                    coverage, percent, qlen, targ, targ, to_keep, to_drop, gene_idx
                )

        # reverse strand
        if current_strand == "strandR":
            if current_scaffold == next_scaffold and current_strand == next_strand:
                ex_targ = targ
                # exon loop
                print('exon_loop')
                while (
                    next_pos_end < current_pos_st + (qlen // 2)
                    and current_scaffold == next_scaffold
                    and current_strand == next_strand
                    and ex_targ < (len(th_df) - 1)):
                   # print(f'cov {coverage} and extarg {ex_targ}')
                    ex_targ += 1
                    (
                        current_scaffold,
                        current_strand,
                        current_pos_st,
                        current_pos_end,
                    ) = _attribute(ex_targ, th_df)
                    if ex_targ != (len(th_df) - 1):
                        (next_scaffold, 
                         next_strand, 
                         next_pos_st, 
                         next_pos_end) = _attribute(ex_targ + 1, th_df)
                        coverage += current_pos_end - current_pos_st
                    if ex_targ == (len(th_df) - 1):
                        coverage += current_pos_end - current_pos_st
                    
                _coverage_test(
                    coverage, percent, qlen, targ, ex_targ, to_keep, to_drop, gene_idx
                )
                targ = ex_targ
                print(100*coverage/qlen)
            else:
                _coverage_test(
                    coverage, percent, qlen, targ, targ, to_keep, to_drop, gene_idx
                )
                print(100*coverage/qlen)
        targ += 1

    # last target
    if (len(th_df) - 1) not in to_keep:
        targ_last = len(th_df) - 1
        coverage_last = th_df["HMM_to"][targ_last] - th_df["HMM_from"][targ_last]
        _coverage_test(
            coverage_last,
            percent,
            qlen,
            targ_last,
            targ_last,
            to_keep,
            to_drop,
            gene_idx,
        )

    # short hits removal
    to_drop = list(dict.fromkeys(to_drop))

    for i in to_drop:
        th_df = th_df.drop([i])

    th_df.reset_index(inplace=True) #Need to keep th_df

    return th_df, gene_idx


def df_to_mat(target_hits_df: pd.DataFrame) -> np.array:
    """
    Produce a matrix of 1 and 0 from a given dataframe.
    """
    qlen = target_hits_df["qlen"][0]
    targ_mat = np.zeros((len(target_hits_df), qlen))
    for tgt in range(len(target_hits_df)):
        st = target_hits_df["HMM_from"][tgt] - 1
        en = target_hits_df["HMM_to"][tgt] - 1
        targ_mat[tgt, st:en] += 1
    return targ_mat


def _parse_sp_fam_name(ofile_path: str) -> tuple[str, str]:
    """
    Parse the species and gene family name for a given file.
    """
    first_name = ofile_path.strip().split("/")[-1]
    names = first_name.strip().split(".")
    # sp_name, fam_name
    return names[1], names[0]


def _raw_count(file, percent_ini, fam_name):
    """
    Initial steps of the count until graph cleaner.
    """
    sorted_df = sort_domtout(file)
    print('sorted')
    print(sorted_df[["target_name", "qlen", "HMM_from", "HMM_to"]])
    target_hits_df = sort_hits(sorted_df)
    print('targets')
    print(target_hits_df[["target_name", "qlen", "HMM_from", "HMM_to"]])
    targetcleaned = graph_cleaner(target_hits_df, percent_ini / 100)
    print('cleaned')
    print(targetcleaned)
    return targetcleaned


def save_clean(targetcleaned, sp_name, fam_name, out):
    try:
        os.mkdir(f'{out}/{sp_name}')
    except OSError as error:
        print(error)
    targetcleaned[0].to_csv(f'{out}/{sp_name}/{fam_name}.{sp_name}.targetcleaned.csv',  index=False, sep=";")
    #print(targetcleaned[1])
    #print(str(targetcleaned[1]))
    with open(f'{out}/{sp_name}/{fam_name}.{sp_name}.geneid.txt', 'w') as gi:
        gi.write(str(targetcleaned[1]))

def cleaning(
    folder: str,
    percent: int,
    outpath: str
) -> str:

    file_list = [
        os.path.join(folder, cfile) for cfile in os.listdir(folder)
    ]

    for file in file_list:
        sp_name, fam_name = _parse_sp_fam_name(file)
        print(file)

        if _file_not_empty(file):
            targetcleaned = _raw_count(file, percent, fam_name)
            save_clean(targetcleaned, sp_name, fam_name, outpath)

    return 1



def main():
    """
    The main function of the program.
    """
    parser = argparse.ArgumentParser(
        description="Clean the domout table and save the cleaned dataframe "
        "and the gene id list in the output folder."
    )
    parser.add_argument(
        "-d",
        "--domtfolder",
        help="The folder in which are the domtout outputs of the hmmsearch.",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-c", 
        "--percent",
        help="The percentage of coverage.",
        type=float,
        default=70
    )

    parser.add_argument(
        "-o",
        "--outfol",
        help="The folder name to the targetcleaned dataframe and list output.",
        type=str,
        default="output/targetcleaned",
    )

    
    args = parser.parse_args()

    try:
        os.mkdir("output/targetcleaned")
    except OSError as error:
        print(error)

    cleaning(
        args.domtfolder,
        args.percent,
        args.outfol,
    )


if __name__ == "__main__":
    main()
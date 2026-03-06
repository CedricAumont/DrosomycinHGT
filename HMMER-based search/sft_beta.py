# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:01:43 2023

@author: caumont
"""
import argparse
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

def strand_choice(seq_file: str) -> str:
    return(f'{seq_file.strip().split("_")[7][-1]}')


def find_orfs_with_trans(seq_file, min_protein_length, strand_to_trans, path_to_seqfile, trans_table=1):
    seq_rec = SeqIO.parse(f'{path_to_seqfile}/{seq_file}', "fasta")
    answer = []
    if strand_to_trans == "A":
        strand_to_trans = strand_choice(seq_file)
    for scaffold, seq_obj in enumerate(seq_rec):
        seq_len = len(seq_obj)
        print(seq_len)
        if strand_to_trans == "B":
                strand_seqobj = (("F", seq_obj.seq),("R", seq_obj.seq.reverse_complement()))
        if strand_to_trans == "F":
                strand_seqobj = (("F", seq_obj.seq),)
        if strand_to_trans == "R":
                strand_seqobj = (("R", seq_obj.seq.reverse_complement()),)
        for strand, nuc in strand_seqobj :
            for frame in range(3):
                trans = nuc[frame:].translate(trans_table)
                trans_len = len(trans)
                aa_start, aa_end = 0, 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end - aa_start >= min_protein_length:
                        if aa_end - aa_start > 100000:
                            print("Too many X", aa_end - aa_start,
                                  trans[aa_start:aa_start+200])
                            aa_end = aa_start + 90000
                        if strand == "F":
                            start = frame + aa_start * 3
                            end = min(seq_len, frame + aa_end * 3 + 3)
                        else:
                            start = seq_len - frame - aa_end * 3 - 3
                            end = seq_len - frame - aa_start * 3
                        answer.append((start, end, strand,
                                       trans[aa_start:aa_end], scaffold+1,
                                       seq_obj.description))
                    aa_start = aa_end + 1
    return answer


def write_results(path_name, seq_file, orf_list):
    gen = f'{seq_file.strip().split("_")[0]}'
    spe = f'{seq_file.strip().split("_")[1]}'
    outfile = f'{path_name}/{gen}_{spe}_protein.faa'
    with open(outfile, 'w') as of:
        for start, end, strand, pro, scaff, descr in orf_list:
            of.write(f'>gac{scaff}_{gen}_{spe}_{start}:{end}, length {len(pro)}, '
                     f'strand{strand}, {start}:{end}, scaffold {scaff}, '
                     f'{descr}')
            of.write(f'\n{pro}\n')


def main():
    """
    The main function of the program.
    """
    parser = argparse.ArgumentParser(description='Six frame translate the target genomes.')
    parser.add_argument("-f", "--seqfile", help="The name of the sequence file to translate.", type=str, required=True)
    parser.add_argument("-p", "--prolen", help="The minimal length for a protein to be translated, default is 15 aa.", type=float, default=15)
    parser.add_argument("-o", "--outpath", help="The path to store results in. Default is protein, but single folder is advised", type=str, default="input/protein/")
    parser.add_argument("-i", "--inpath", help="The path to seqfile.", type=str, default=".")
    parser.add_argument("-s", "--strand", help="The strand translated, F: +1, R: -1, B: Both, A: Auto", type=str, default="B")
    args = parser.parse_args()

    orf_list = find_orfs_with_trans(args.seqfile, args.prolen, args.strand, args.inpath)
    write_results(args.outpath, args.seqfile, orf_list)


if __name__ == "__main__":
    main()
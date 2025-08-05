# -*- coding: utf-8 -*-
"""
@created: Jun 02 2023
@modified: Jul 02 2024
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Apply additional filtering on the tables of all mutations.
"""

import argparse
import os
import gzip
import numpy  as     np
import pandas as     pd
import re
import subprocess

# functions ============================================================================================================

def read_header(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as file:
            header = [x for x in file.readlines() if x.startswith("##")]
    else:
        with open(path, "r") as file:
            header = [x for x in file.readlines() if x.startswith("##")]
    return header


def read_table(path, **kwargs):
    header = read_header(path)
    df = pd.read_table(path, skiprows=len(header), na_values=["-","."], **kwargs)
    return df


def save_table_with_header(df, header, output):
    output_header = output + "header.tsv"

    print("-writing header in file %s" % output_header)
    with open(output_header, "w") as file_header:
        for line in header:
            file_header.write(line)

    # write contents
    if output.endswith(".gz"):
        output_uncompressed = output[:-3] + ".tmp"
        output_concatenate = output[:-3]
    else:
        output_uncompressed = output + ".tmp"
        output_concatenate = output

    print("-writing contents in file %s" % output_uncompressed)
    df.to_csv(output_uncompressed, index=False, sep="\t")

    # # concat both files
    cmd = "cat %s %s >> %s" % (output_header, output_uncompressed, output_concatenate)
    print("-running the command:\n\t%s" % cmd)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    os.remove(output_header)
    os.remove(output_uncompressed)

    # compress if required
    if output_concatenate != output:
        cmd = "gzip %s" % output_concatenate
        print("-running the command:\n\t%s" % cmd)
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def filter_regex(df, regex, column, mode="out", return_mask=False):
    mask_regex = df[column].str.contains(regex)
    mask_regex = mask_regex.replace(np.nan, False)
    if mode=="out":
        print("-filtered out %d/%d mutations using column %s and regex %s" % \
              (sum(mask_regex), mask_regex.shape[0], column, regex))
    elif mode=="in":
        mask_regex = ~mask_regex
        print("-filtered in %d/%d mutations using column %s and regex %s" % \
              (sum(~mask_regex), mask_regex.shape[0], column, regex))
    else:
        raise ValueError("-ERROR: unrecognized value of mode %s" % mode)

    if return_mask:
        return ~mask_regex
    else:
        return df.loc[~mask_regex].copy()


def add_filter_not_exonic(df):
    col_filt = "FILTER"
    col = "Variant_Classification"
    filt = "not_exonic"
    vals_keep = ["3'UTR", "5'UTR", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent",
                 "Splice_Site", "Translation_Start_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                 "In_Frame_Ins"]
    regex_keep = "|".join(vals_keep)

    mask_keep = filter_regex(df, regex_keep, col, mode="in", return_mask=True)
    mask_pass = df[col_filt]=="PASS"
    df.loc[~mask_keep & mask_pass, col_filt] = filt
    df.loc[~mask_keep & ~mask_pass, col_filt] += ",%s" % filt

    return df


def main(args):
    # load data
    df_inp_maf = read_table(args.inp_maf, low_memory=False)
    df_inp_maf_ann = read_table(args.inp_maf_ann, low_memory=False)
    df_inp_maf_civ = read_table(args.inp_maf_civ, low_memory=False)
    df_inp_maf_okb = read_table(args.inp_maf_okb, low_memory=False)

    # refiltering
    print("INFO: adding filter not exonic on maf table...")
    df_inp_maf = add_filter_not_exonic(df_inp_maf)

    print("INFO: adding filter not exonic on maf ann table...")
    df_inp_maf_ann = add_filter_not_exonic(df_inp_maf_ann)

    print("INFO: adding filter not exonic on maf civ table...")
    df_inp_maf_civ = add_filter_not_exonic(df_inp_maf_civ)

    print("INFO: adding filter not exonic on maf okb table...")
    df_inp_maf_okb = add_filter_not_exonic(df_inp_maf_okb)

    # save refiltered tables
    print("INFO: saving df_inp_maf to %s ..." % args.out_maf)
    header_inp_maf = read_header(args.inp_maf)
    save_table_with_header(df_inp_maf, header_inp_maf, args.out_maf)

    print("INFO: saving df_inp_maf_ann to %s ..." % args.out_maf_ann)
    df_inp_maf_ann.to_csv(args.out_maf_ann, sep="\t", index=False)

    print("INFO: saving df_inp_maf_civ to %s ..." % args.out_maf_civ)
    df_inp_maf_civ.to_csv(args.out_maf_civ, sep="\t", index=False)

    print("INFO: saving df_inp_maf_okb to %s ..." % args.out_maf_okb)
    df_inp_maf_okb.to_csv(args.out_maf_okb, sep="\t", index=False)


if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Apply last mutation filters.')
    parser.add_argument('--inp_maf', type=str, help='Path to all currently PASS mutations.')
    parser.add_argument('--inp_maf_ann', type=str, help='Path to all oncogenic mutations.')
    parser.add_argument('--inp_maf_civ', type=str, help='Path to all CIViC-annotated mutations.')
    parser.add_argument('--inp_maf_okb', type=str, help='Path to all OncoKB-annotated mutations.')
    parser.add_argument('--out_maf', type=str, help='Path to all refiltered mutations.')
    parser.add_argument('--out_maf_ann', type=str, help='Path to all refiltered and oncogenic mutations.')
    parser.add_argument('--out_maf_civ', type=str, help='Path to all refiltered and CIViC-annotated mutations.')
    parser.add_argument('--out_maf_okb', type=str, help='Path to all refiltered and OncoKB-annotated mutations.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")
    main(args)

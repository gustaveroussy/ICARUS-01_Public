# -*- coding: utf-8 -*-
"""
@created: Jun 01 2023
@modified: Jun 01 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Convert a VCF into a TSV file of positions dealing with mutliallelic sites.
"""

import argparse
import os
import gzip
import pandas as     pd

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
    df = pd.read_table(path, skiprows=len(header), na_values=None, **kwargs)
    return df


def convert_num_to_str(x, digits=None):
    try:
        if int(x)==x:
            y = "%d" % x
        else:
            raise
    except:
        try:
            if digits is not None:
                y = ("%%.%df" % digits) % x
            else:
                y = str(x)
            if y=="nan":
                y = x
        except:
            y = x
    return y


def shift_multiallelic_pos(x):
    pos = x["POS"]
    ref = x["REF"]
    alt = x["ALT"]
    while len(ref)>0 and len(alt)>0 and ref[0]==alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    if ref=="" or alt=="":
        pos = pos - 1
    shift = pos-x["POS"]
    if shift!=0:
        print("shifted multiallelic site %s at %s:%s by %d" % (x["ALT_Original"], x["#CHROM"], x["POS"], shift))
    return pos


def main(args):
    # load data
    df_inp = read_table(args.inp, low_memory=False)
    print("INFO: loaded %d sites" % df_inp.shape[0])

    # select columns
    cols = ["#CHROM", "POS", "REF", "ALT"]
    df_inp = df_inp[cols].copy()

    # deal with multiallelic sites
    df_inp["ALT_Original"] = df_inp["ALT"]
    mask_multi = df_inp["ALT"].str.contains(",")
    print("INFO: detected %d/%d sites that are multi-allelic" % (sum(mask_multi), len(mask_multi)))

    # select the first alternative allele
    df_inp.loc[mask_multi, "ALT"] = df_inp.loc[mask_multi, "ALT"].apply(lambda x: x.split(",")[0])

    # trim the leading base and shift the position by +1 until REF and ALT do not share the same leading base 
    # if the REF becomes empty, we had an insertion and undo one shift
    # if the ALT becomes empty, we had a deletion and undo one shift
    # otherwise, we have a SNV or MNV
    df_inp.loc[mask_multi, "POS"] = df_inp.loc[mask_multi, :].apply(shift_multiallelic_pos, axis=1)

    # save
    df_out = df_inp[["#CHROM", "POS"]].drop_duplicates()
    for col in df_out:
        df_out[col] = df_out[col].apply(convert_num_to_str)
    print("INFO: saving %d positions" % df_out.shape[0])
    df_out.to_csv(args.out, sep="\t", index=False, header=False)


if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Convert VCF into a TSV file of positions.')
    parser.add_argument('--inp', type=str, help='Path to input.')
    parser.add_argument('--out', type=str, help='Path to output.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")
    main(args)

# -*- coding: utf-8 -*-
"""
@created: May 13 2023
@modified: Jun 02 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Add pileup data of the tumor sample on the annotated positions from the normal sample.
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
    df = pd.read_table(path, skiprows=len(header), na_values=None, **kwargs)
    return df


def save_table_with_header(df, header, output):
    output_header = output + ".header.tsv"

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


def get_depth(x):
    try:
        if x=="*":
            return 1
        else:
            # remove ^ and following char indicating read end and mapping quality
            y = re.sub("(?<=\\^)(.)", "", x)
            y = re.sub("\\^", "" , y)

            # remove $ char indicating read end
            y = re.sub("\\$", "" , y)

            # count ins sizes
            ins_sizes = re.findall("\+[0-9]+", y)
            ins_sizes = sum([abs(int(s))+len(s) for s in ins_sizes])

            # count del sizes
            del_sizes = re.findall("-[0-9]+", y)
            del_sizes = sum([abs(int(s))+len(s) for s in del_sizes])

            return len(list(y)) - ins_sizes - del_sizes

    except:
        print(x)
        return 0


def get_ref_count(x):
    try:
        # remove ^ and following char indicating read end and mapping quality
        y = re.sub("(?<=\\^)(.)", "", x)
        y = re.sub("\\^", "" , y)

        # remove $ char indicating read end
        y = re.sub("\\$", "" , y)

        # count ins sizes
        ins_sizes = re.findall("\+[0-9]+", y)

        # count del sizes
        del_sizes = re.findall("-[0-9]+", y)

        return y.count(".") + y.count(",") - len(del_sizes) - len(ins_sizes)
    except:
        print(x)
        return 0


def process_ins(read_bases):
    # extract insertions
    ins_sizes = re.findall("\+[0-9]+", read_bases)
    ins_sizes = [abs(int(s)) for s in ins_sizes]
    read_bases_no_ins = read_bases
    ins_all = []
    ins_pos = []
    for ins_size in ins_sizes:
        regex_ins = "\+%d[ACGTNacgtn*#]{%d}" % (ins_size, ins_size)
        ins_one = re.findall(regex_ins, read_bases_no_ins)[0]
        ins_one = ins_one.replace("+%d" % ins_size, "")
        pos = re.search(regex_ins, read_bases_no_ins).start()-1
        read_bases_no_ins = re.sub(regex_ins, "", read_bases_no_ins, count=1)
        ins_all.append(ins_one)
        ins_pos.append(pos)
    return read_bases_no_ins, ins_all, ins_pos


def process_del(read_bases):
    # extract deletions
    del_sizes = re.findall("-[0-9]+", read_bases)
    del_sizes = [abs(int(s)) for s in del_sizes]
    read_bases_no_del = read_bases
    del_all = []
    del_pos = []
    for del_size in del_sizes:
        regex_del = "-%d[ACGTNacgtn*#]{%d}" % (del_size, del_size)
        del_one = re.findall(regex_del, read_bases_no_del)[0]
        del_one = del_one.replace("-%d" % del_size, "")
        pos = re.search(regex_del, read_bases_no_del).start()-1
        read_bases_no_del = re.sub(regex_del, "", read_bases_no_del, count=1)
        del_all.append(del_one)
        del_pos.append(pos)
    return read_bases_no_del, del_all, del_pos


def get_alt_count(x, base_qual=False, mapping_qual=False):
    mut_type =  x["Variant_Type"]
    read_bases = x["t2_read_bases"]
    if base_qual:
        qual_bases = x["t2_qual_bases"]
    if mapping_qual:
        qual_mapping = x["t2_qual_mapping"]

    if type(read_bases)!=str:
        return np.nan
    else:
        # remove ^ and following char indicating read end and mapping quality
        read_bases = re.sub("(?<=\\^)(.)", "", read_bases)
        read_bases = re.sub("\\^", "" , read_bases)

        # remove $ char indicating read end
        read_bases = re.sub("\\$", "" , read_bases)

        # order of processing of insertion/deletion is important to get correct indices for qualities
        if mut_type=="DEL" or mut_type!="INS":
            read_bases_no_ins, ins_all, ins_all_pos = process_ins(read_bases)
            read_bases_clean, del_all, del_pos = process_del(read_bases_no_ins)
        elif mut_type=="INS":
            read_bases_no_del, del_all, del_all_pos = process_del(read_bases)
            read_bases_clean, ins_all, ins_pos = process_ins(read_bases_no_del)

        if mut_type not in["DEL", "INS"]:
            alt = x["Tumor_Seq_Allele2"]
            count = read_bases_clean.upper().count(alt)
            alt_pos = [m.start(0) for m in re.finditer(alt, read_bases_clean.upper())]
        elif mut_type=="INS":
            alt = x["Tumor_Seq_Allele2"]
            count = len([e for e in ins_all if e.upper()==alt])
            alt_pos = [i for i,e in zip(ins_pos, ins_all) if e.upper()==alt]
        elif mut_type=="DEL":
            ref = x["Reference_Allele"]
            count = len([e for e in del_all if e.upper()==ref])
            alt_pos = [i for i,e in zip(del_pos, del_all) if e.upper()==ref]

    if base_qual and qual_bases is not None and read_bases!="*":
        qual_b = ";".join([str(ord(qual_bases[pos])-33) for pos in alt_pos])
    else:
        qual_b = None
    if mapping_qual and qual_mapping is not None and read_bases!="*":
        qual_m = ";".join([str(ord(qual_mapping[pos])-33) for pos in alt_pos])
    else:
        qual_m = None

    if base_qual and mapping_qual:
        return count, qual_b, qual_m
    elif base_qual:
        return count, qual_b
    elif mapping_qual:
        return count, qual_m
    else:
        return count


def load_and_process_mpileup(path):
    df_pup = read_table(path, header=None, keep_default_na=False, quoting=3)
    cols = ["Chromosome", "Start_Position", "Reference", "t2_depth", "t2_read_bases", "t2_qual_bases"]
    if df_pup.shape[1] == len(cols)+1:
        cols += ["Tumor_Sample_Barcode_2"]
    elif df_pup.shape[1] == len(cols)+2:
        cols += ["t2_qual_mapping", "Tumor_Sample_Barcode_2"]
    df_pup.columns = cols

    # add t2_ref_count
    assert df_pup["t2_depth"].equals(df_pup["t2_read_bases"].apply(get_depth))
    df_pup["t2_ref_count"] = df_pup["t2_read_bases"].apply(get_ref_count)
    assert all(df_pup["t2_ref_count"] <= df_pup["t2_depth"])

    return df_pup


def main(args):
    # load data
    df_inp_maf = read_table(args.inp_maf, low_memory=False)
    df_inp_pup = load_and_process_mpileup(args.inp_pup)

    # for maf deletions with empty alternative allele, shift position by -1
    mask_del = (df_inp_maf["Variant_Type"]=="DEL") & (df_inp_maf["Tumor_Seq_Allele2"].isnull())
    df_inp_maf.loc[mask_del, "Start_Position"] -= 1

    # cols order
    cols_maf_order = df_inp_maf.columns.tolist()

    # add Order to preserver original order
    df_inp_maf["Order"] = np.arange(df_inp_maf.shape[0])

    # merge on chr/position
    col_row_id = "Row_Id"
    cols_maf_id = ["Chromosome", "Start_Position"]
    cols_pup_id = ["Chromosome", "Start_Position"]

    df_inp_maf[col_row_id] = df_inp_maf[cols_maf_id].astype(str).apply("/".join, axis=1)
    df_inp_pup[col_row_id] = df_inp_pup[cols_pup_id].astype(str).apply("/".join, axis=1)

    assert df_inp_pup.shape[0]==df_inp_pup[col_row_id].nunique()

    cols_inp_pup = [col_row_id, "Tumor_Sample_Barcode_2"] + [x for x in df_inp_pup if x.startswith("t2_")]
    n_row_inp = df_inp_maf.shape[0]
    df_out_maf = df_inp_maf.merge(df_inp_pup[cols_inp_pup], how="left", on=col_row_id)
    n_row_out = df_out_maf.shape[0]
    assert n_row_inp==n_row_out

    # add t_alt_count
    index_t_alt_count = df_out_maf.columns.get_loc("t2_ref_count") + 1
    df_count = df_out_maf.apply(get_alt_count, axis=1, base_qual=False)
    df_count = df_count.apply(pd.Series)
    df_count.columns = ["t2_alt_count"]
    df_out_maf.insert(index_t_alt_count, "t2_alt_count", df_count["t2_alt_count"])

    # reset position for deletions 
    mask_del = (df_out_maf["Variant_Type"]=="DEL") & (df_out_maf["Tumor_Seq_Allele2"].isnull())
    df_out_maf.loc[mask_del, "Start_Position"] += 1

    # count missed positions
    ids_pup_not_maf = list(set(df_inp_pup[col_row_id]).difference(set(df_inp_maf[col_row_id])))
    if len(ids_pup_not_maf):
        print("-WARNING: %d positions in pileup do not match any position in the MAF table" % len(ids_pup_not_maf))
        if len(ids_pup_not_maf)<10:
            print("\t" + "\n\t".join(ids_pup_not_maf))

    ids_maf_not_pup = list(set(df_inp_maf[col_row_id]).difference(set(df_inp_pup[col_row_id])))
    if len(ids_maf_not_pup):
        print("-WARNING: %d positions in MAF table do not match any position in the pileup" % len(ids_maf_not_pup))
        if len(ids_maf_not_pup)<10:
            print("\t" + "\n\t".join(ids_maf_not_pup))

    # format
    df_out_maf = df_out_maf.sort_values(by="Order")
    del df_out_maf["Order"]
    del df_out_maf[col_row_id]
    df_out_maf["HGVSp"] = df_out_maf["HGVSp"].str.replace("%3D", "=")
    df_out_maf["HGVSp_Short"] = df_out_maf["HGVSp_Short"].str.replace("%3D", "=")
    df_out_maf["all_effects"] = df_out_maf["all_effects"].str.replace("%3D", "=")

    for col in df_out_maf:
        df_out_maf[col] = df_out_maf[col].apply(convert_num_to_str)

    # select columns
    cols = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type",
            "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Tumor_Sample_Barcode_2", "t_depth", "t_alt_count",
            "t_ref_count", "n_depth", "n_alt_count", "n_ref_count"] + [x for x in df_out_maf if x.startswith("t2_")]
    cols = [x for x in cols if x in df_out_maf]

    # save
    print("INFO: saving df_out_maf to %s ..." % args.out_maf)
    df_out_maf[cols].to_csv(args.out_maf, sep="\t", index=False)



if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Apply pileup data to mutations in MAF format.')
    parser.add_argument('--inp_maf', type=str, help='Path to annotated mutations of blood sample.')
    parser.add_argument('--inp_pup', type=str, help='Path to mpileup data of tumor sample.')
    parser.add_argument('--out_maf', type=str, help='Path to annotated mutations of blood sample with tumor pileup.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")
    main(args)

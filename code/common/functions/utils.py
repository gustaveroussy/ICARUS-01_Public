# -*- coding: utf-8 -*-
"""
@created: Jul 07 2022
@modified: Dec 06 2024
@author: Yoann Pradat

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Useful functions for preparing tables for oncoplot-like figures.
"""

import gzip
import numpy as np
import os
import pandas as pd
import re


def read_excel_with_dates(path, skiprows=None, header=0, dtype=None, cols_date=None, dayfirst=True, format=None, **kwargs):
    if dtype is None and cols_date is not None:
        dtype = {k: "str" for k in cols_date}
    elif dtype is not None and cols_date is not None:
        for col in cols_date:
            dtype[col] = "str"
    elif cols_date is not None:
        raise ValueError("If using both cols_date and dtype arguments, please specify a dictionary for dtype.")

    df_cln = pd.read_excel(path, na_values=["na", "na ", "NA", ".", "not available"], header=header,
                           skiprows=skiprows, parse_dates=False, engine="openpyxl", dtype=dtype, **kwargs)

    if cols_date is not None:
        for col in cols_date:
            if format is not None:
                df_cln.loc[:, col] = pd.to_datetime(df_cln[col], format=format, errors="coerce")
            else:
                df_cln.loc[:, col] = pd.to_datetime(df_cln[col], dayfirst=dayfirst, errors="coerce")

    return df_cln


def read_header(path, prefix):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    else:
        with open(path, "r") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    return header


def load_table(path, header=0, dtype=None, header_prefix=None, **kwargs):
    if header_prefix is not None:
        skiprows = len(read_header(path=path, prefix=header_prefix))
    else:
        skiprows = None

    if path.endswith("xlsx"):
        df = read_excel_with_dates(path, skiprows=skiprows, header=header, dtype=dtype, **kwargs)
    elif ".txt" in path or ".tsv" in path or ".maf" in path:
        if "sep" not in kwargs:
            df = pd.read_csv(path, skiprows=skiprows, header=header, dtype=dtype, sep="\t", **kwargs)
        else:
            df = pd.read_csv(path, skiprows=skiprows, header=header, dtype=dtype, **kwargs)
    else:
        df = pd.read_csv(path, skiprows=skiprows, header=header, dtype=dtype, **kwargs)

    return df


def join_unique(x):
    return "|".join(list(set(x)))


def concatenate_cna(df_alt, df_cna, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim):
    if df_cna.shape[0]==0:
        return df_alt
    else:
        cols_ids = ["Aliquot_Id", "Sample_Id", "Cluster_Id", "Subject_Id"]
        cols_ids = [x for x in cols_ids if x in df_cna]

        cols_det = ["Copy_Number_More", "TCN_EM:LCN_EM"]
        cols = cols_ids + [col_gen, col_alt_cat, col_alt_cat_sim, col_alt, col_alt_det, "Cellular_Fraction", "Annotated"]
        df_cna["Cellular_Fraction"] = df_cna["cf.em"]
        df_cna["Annotated"] = "No"
        df_cna[col_gen] = df_cna["Hugo_Symbol"]
        mask_del = df_cna["Copy_Number"]==-2
        df_cna.loc[mask_del, col_alt_cat] = "Deletion"
        df_cna.loc[mask_del, col_alt_cat_sim] = "Del"
        df_cna.loc[mask_del, col_alt] = "Del"
        df_cna.loc[mask_del, col_alt_det] = df_cna.loc[mask_del, cols_det].apply(" - ".join, axis=1)

        mask_amp = df_cna["Copy_Number"]==2
        df_cna.loc[mask_amp, col_alt_cat] = "Amplification"
        df_cna.loc[mask_amp, col_alt_cat_sim] = "Amp"
        df_cna.loc[mask_amp, col_alt] = "Amp"
        df_cna.loc[mask_amp, col_alt_det] = df_cna.loc[mask_amp, cols_det].apply(" - ".join, axis=1)

        mask_oth = (~mask_del & ~mask_amp) & (~df_cna["Copy_Number_More"].isnull())
        df_cna.loc[mask_oth, col_alt_det] = df_cna.loc[mask_oth, cols_det].apply(" - ".join, axis=1)
        df_cna.loc[mask_oth, col_alt] = df_cna.loc[mask_oth, col_alt_det]

        mask_null = df_cna[col_alt_cat].isnull()
        df_cna.loc[mask_null, col_alt_cat] = df_cna.loc[mask_null, "Copy_Number_More"]

        df_alt = pd.concat((df_alt, df_cna[cols]), axis=0)
        df_alt = df_alt.drop_duplicates(subset=[col_sam_id, col_gen, col_alt_det], keep="first")

        return df_alt


def concatenate_mut(df_alt, df_mut, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim,
                    mut_det_level):
    if df_mut.shape[0]==0:
        return df_alt
    else:
        cols_ids = ["Aliquot_Id", "Sample_Id", "Cluster_Id", "Subject_Id"]
        cols_ids = [x for x in cols_ids if x in df_mut]
        cols = cols_ids + [col_gen, "Variant_Classification", col_alt_cat, col_alt_cat_sim, col_alt, col_alt_det,
                           "Annotated", "t_vaf"]

        df_mut["Annotated"] = "No"
        df_mut[col_gen] = df_mut["Hugo_Symbol"]
        if mut_det_level=="HGVSp":
            get_mut_det = lambda x: x["HGVSp_Short"] if type(x["HGVSp_Short"])==str else x["HGVSc"]
        elif mut_det_level=="HGVSc":
            get_mut_det = lambda x: x["HGVSc"]
        else:
            raise ValueError("Unsupported value %s for 'mut_det_level'. Choose 'HGVSp' or 'HGVSc'" % mut_det_level)
        hgvsp_short_split = lambda x: x.split("p.")[1].replace("%3D", "=") if type(x)==str else x
        exon_split = lambda x: x.split("/")[0] if type(x)==str else ""

        mask_ins = df_mut["Variant_Classification"].isin(["Frame_Shift_Ins", "In_Frame_Ins"])
        df_mut.loc[mask_ins, col_alt_cat] = "Indel"
        df_mut.loc[mask_ins, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_ins, col_alt] = "Exon " + df_mut.loc[mask_ins,"EXON"].apply(exon_split) + " Ins"
        df_mut.loc[mask_ins, col_alt_det] = df_mut.loc[mask_ins, ["HGVSp_Short", "HGVSc"]].apply(get_mut_det, axis=1)

        mask_del = df_mut["Variant_Classification"].isin(["Frame_Shift_Del", "In_Frame_Del"])
        df_mut.loc[mask_del, col_alt_cat] = "Indel"
        df_mut.loc[mask_del, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_del, col_alt] = "Exon " + df_mut.loc[mask_del,"EXON"].apply(exon_split) + " Del"
        df_mut.loc[mask_del, col_alt_det] = df_mut.loc[mask_del, ["HGVSp_Short", "HGVSc"]].apply(get_mut_det, axis=1)

        mask_mut = (~mask_ins) & (~mask_del)
        mask_nul_hgvsp = df_mut["HGVSp_Short"].isnull()
        df_mut.loc[mask_mut, col_alt_cat] = "Mutation"
        df_mut.loc[mask_mut, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_mut, col_alt] = df_mut.loc[mask_mut, "HGVSp_Short"].apply(hgvsp_short_split)
        df_mut.loc[mask_mut & mask_nul_hgvsp, col_alt] = \
                df_mut.loc[mask_mut & mask_nul_hgvsp, "Variant_Classification"]
        df_mut.loc[mask_mut, col_alt_det] = df_mut.loc[mask_mut, ["HGVSp_Short", "HGVSc"]].apply(get_mut_det, axis=1)

        # for splice, set col_alt to Splice_Site
        mask_splice = df_mut["Variant_Classification"]=="Splice_Site"
        df_mut.loc[mask_splice, col_alt] = "Splice_Site"
        mask_splice = df_mut[col_alt_det].apply(lambda x: "splice" in x if type(x)==str else False)
        df_mut.loc[mask_splice, col_alt_det] = "Splice_Site"

        # fill null col_alt_det
        mask_null = df_mut[col_alt_det].isnull()

        if mut_det_level=="HGVSp":
            df_mut.loc[mask_null, col_alt_det] = df_mut.loc[mask_null, "Variant_Classification"]
        elif mut_det_level=="HGVSc":
            df_mut["Mutation_Id"] = df_mut["Start_Position"].astype(str) + ":" + \
                df_mut[["Reference_Allele", "Tumor_Seq_Allele2"]].fillna("").apply(">".join, axis=1)
            df_mut.loc[mask_null, col_alt_det] = df_mut.loc[mask_null, "Mutation_Id"]

        df_alt = pd.concat((df_alt, df_mut[cols]), axis=0)
        df_alt = df_alt.drop_duplicates(subset=[col_sam_id, col_gen, col_alt_det], keep="first")

        return df_alt


def process_alt(df_alt, df_cna, df_mut, sam_list, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                keep_alt_det=False, keep_alt_det_one=True, mut_det_level="HGVSp"):

    # create col_alt_cat_sim for naming the alerations and group some categories for coloring
    col_alt_cat_sim = "%s_Simple" % col_alt_cat

    if df_alt is not None and len(df_alt)>0:
        df_alt[col_alt_cat_sim] = df_alt[col_alt_cat].replace({"Del": "Mut", "Ins": "Mut", "Amplification": "Amp",
                                                               "Deletion": "Del"})
        df_alt[col_alt_cat] = df_alt[col_alt_cat].replace({"Mut": "Mutation", "Ins": "Indel", "Del": "Indel"})
        df_alt["Annotated"] = "Yes"
        df_alt = df_alt.loc[df_alt[col_sam_id].isin(sam_list)].copy()
        df_alt = df_alt.rename(columns={"Tumor_Type": "Cancer_Type"})
    if df_cna is not None:
        df_cna = df_cna.loc[df_cna[col_sam_id].isin(sam_list)].copy()
    if df_mut is not None:
        df_mut = df_mut.loc[df_mut[col_sam_id].isin(sam_list)].copy()

    # concatenate cna
    if df_cna is not None:
        df_alt = concatenate_cna(df_alt, df_cna, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                                 col_alt_cat_sim)
        df_alt = df_alt.reset_index(drop=True)

    # concatenate mut
    if df_mut is not None:
        df_alt = concatenate_mut(df_alt, df_mut, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                                 col_alt_cat_sim, mut_det_level)
        df_alt = df_alt.reset_index(drop=True)

    # there should not be any null col_alt_det
    assert df_alt[col_alt_det].isnull().sum()==0

    # concatenate gene name and alteration
    # for add gene symbol as prefix except for
    #  - alterations that are already equal to gene symbol
    df_gen_alt = df_alt[[col_gen, col_alt]].drop_duplicates().fillna("NA")
    df_gen_alt_gby = df_gen_alt.groupby(col_gen)[col_alt].agg("|".join)
    mask_a = (df_alt[col_alt] != df_alt[col_gen])
    df_alt.loc[mask_a, col_alt] = df_alt.loc[mask_a, [col_gen, col_alt]].apply(lambda x: " ".join(x.dropna()), axis=1)

    # if don't keep alteration detail, only the gene name is kept unless only 1 alt is associated with the gene,
    # in which case replace by alt
    if not keep_alt_det:
        df_gen_alt = df_alt[[col_gen, col_alt]].drop_duplicates()
        df_gen_alt_gby = df_gen_alt.groupby(col_gen).size()
        if keep_alt_det_one:
            gen_1_alt_only = df_gen_alt_gby[df_gen_alt_gby==1].index.tolist()
            mask_1_alt_only = df_alt[col_gen].isin(gen_1_alt_only)
            df_alt.loc[~mask_1_alt_only,col_alt] = df_alt.loc[~mask_1_alt_only, col_gen]
        else:
            df_alt[col_alt] = df_alt[col_gen]

    return df_alt


def remove_flagged_mutations(df_mut, df_mut_flag):
    cols_mut = df_mut.columns.tolist()
    col_chr = "Chromosome"
    col_pos = "Start_Position"
    col_ref, col_alt = "Reference_Allele", "Tumor_Seq_Allele2"
    cols_id = [col_chr, col_pos, col_ref, col_alt]
    col_id = "Mutation_Id"
    df_mut[col_id] = df_mut[cols_id].fillna("N/A").astype(str).apply("/".join, axis=1)
    df_mut_flag[col_id] = df_mut_flag[cols_id].fillna("N/A").astype(str).apply("/".join, axis=1)

    # remove according to decision
    mask_remove_all = df_mut_flag["Decision"] == "REMOVE"
    mask_remove_unm = df_mut_flag["Decision"] == "REMOVE IN UNMATCHED"
    ids_remove_all = df_mut_flag.loc[mask_remove_all, col_id].tolist()
    ids_remove_unm = df_mut_flag.loc[mask_remove_unm, col_id].tolist()
    mask_mut_remove_all = df_mut[col_id].isin(ids_remove_all)
    mask_mut_remove_unm = df_mut[col_id].isin(ids_remove_all) & df_mut["Matched_Norm_Sample_Barcode"].isnull()
    mask_keep = (~mask_mut_remove_all) & (~mask_mut_remove_unm)

    print("-dropped %d/%d flagged mutations" % (sum(~mask_keep), mask_keep.shape[0]))
    return df_mut.loc[mask_keep, cols_mut].copy()


def combine_all_alterations(alt, cln, cna, mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_sub_id, col_sam_id,
                            mut_flag=None, seqs_select="dna", subs_select="all", sams_select="all", keep_alt_det=False,
                            keep_alt_det_one=True, non_syn_mut_only=True, cna_selection="focal_non_neutral",
                            mut_det_level="HGVSp"):

    # get clinical table
    df_cln = pd.read_table(cln)

    # reorganize clinical table per visit
    visits = ["BAS", "ONT", "EOT"]
    mask_b = df_cln["Cancer_Type"] == "Breast Cancer"
    mask_l = df_cln["Cancer_Type"] == "Lung Cancer"
    cols_oth = [x for x in df_cln if not any([v in x for v in visits])]
    dfs_cln = []
    for visit in visits:
        cols_visit = [x for x in df_cln if visit in x]
        if len(cols_visit)>0:
            cols_visit_id = [x for x in cols_visit if x.startswith("Sample_Id")]
            df_cln_visit = df_cln[cols_oth + cols_visit]
            df_cln_visit = df_cln_visit.loc[df_cln_visit[cols_visit_id].isnull().mean(axis=1)<1].copy()
            df_cln_visit = df_cln_visit.rename(columns={x: x.replace("_%s" % visit, "") for x in cols_visit})
            df_cln_visit.loc[mask_b, "Cohort"] = "BREAST-%s" % visit
            df_cln_visit.loc[mask_l, "Cohort"] = "LUNG-%s" % visit
            df_cln_visit["Timepoint"] = visit
            dfs_cln.append(df_cln_visit)
    df_cln = pd.concat(dfs_cln)
    df_cln["Sample_Id_DNA_P"] = df_cln[["Sample_Id_DNA_T", "Sample_Id_DNA_N"]].fillna("NA").apply("_vs_".join, axis=1)
    mask_dna = ~df_cln["Sample_Id_DNA_T"].isnull()
    df_cln.loc[mask_dna, col_sam_id] = df_cln.loc[mask_dna, "Sample_Id_DNA_T"].str[:-4]
    # TODO
    if sum(~mask_dna)>0:
        df_cln.loc[~mask_dna, col_sam_id] = df_cln.loc[~mask_dna, "Sample_Id_RNA_T"].str[:-4]

    # get alterations table if any
    if alt is not None:
        df_alt = pd.read_table(alt)
    else:
        df_alt = None

    if cna is not None:
        df_cna = pd.read_table(cna)
        if cna_selection=="focal_high_level":
            mask_keep = ~df_cna["Copy_Number"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are focal and HL/ML or HD" % (sum(mask_keep), mask_keep.shape[0]))
        elif cna_selection=="focal_non_neutral":
            mask_keep = ~df_cna["Copy_Number_More"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are focal and non-neutral" % (sum(mask_keep), mask_keep.shape[0]))
        elif cna_selection=="focal_all":
            print("-retained all %d focal CNAs (some may be neutral)" % df_cna.shape[0])
        else:
            choices = ["focal_high_level", "focal_non_neutral", "focal_all"]
            raise ValueError("Unknown value '%s' for cna_selection. Choose one of: %s" % (cna_selection, choices))
    else:
        df_cna = None

    if mut is not None:
        df_mut = pd.read_table(mut, sep="\t", skiprows=len(read_header(path=mut, prefix="##")), low_memory=False)
        df_mut["HGVSp_Short"] = df_mut["HGVSp_Short"].apply(lambda x: x.replace("%3D", "=") if type(x)==str else x)
        if non_syn_mut_only:
            mask_keep = df_mut["Variant_Classification"]!="Silent"
            df_mut = df_mut.loc[mask_keep].copy()
            print("-retained %d/%d mutations that non-synonymous" % (sum(mask_keep), mask_keep.shape[0]))
        if mut_flag is not None:
            df_mut_flag = pd.read_excel(mut_flag, dtype=str)
            df_mut = remove_flagged_mutations(df_mut, df_mut_flag)
    else:
        df_mut = None

    # add Sample_Id to all tables
    col_tsb = "Tumor_Sample_Barcode"
    col_nsb = "Matched_Norm_Sample_Barcode"
    cols_cln_dna = ["Sample_Id_DNA_P","Sample_Id","Subject_Id"]
    # TODO
    cols_cln_rna = ["Sample_Id_RNA_T","Sample_Id","Subject_Id"]

    if df_alt is not None:
        df_alt["Sample_Id"] = df_alt["Biopsy_Id"]
        del df_alt["Biopsy_Id"]

    if df_cna is not None:
        df_cna["Sample_Id_DNA_P"] = df_cna[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
        df_cna = df_cna.merge(df_cln[cols_cln_dna], how="left", on="Sample_Id_DNA_P")

    if df_mut is not None:
        df_mut["Sample_Id_DNA_P"] = df_mut[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
        df_mut = df_mut.merge(df_cln[cols_cln_dna], how="left", on="Sample_Id_DNA_P")

        # for mut, add t_vaf
        df_mut["t_vaf"] = df_mut["t_alt_count"]/df_mut["t_depth"]

    # select only samples that pass QC unless specified otherwise by the user
    if "QC_DNA" in df_cln:
        mask_dna = df_cln["QC_DNA"]=="PASS"
    else:
        mask_dna = (~df_cln["Sample_Id_DNA_T"].isnull())

    # TODO
    # if "QC_RNA" in df_cln:
    #     mask_rna = df_cln["QC_RNA"]=="PASS"
    # else:
    #     mask_rna = (~df_cln["Sample_Id_RNA_T"].isnull())

    if seqs_select.lower()=="dna":
        mask_qc = mask_dna
    elif seqs_select.lower()=="rna":
        mask_qc = mask_rna
    elif seqs_select.lower()=="dnarna":
        # select only samples with both dna and rna
        sams_dna = df_cln.loc[mask_dna, col_sam_id]
        sams_rna = df_cln.loc[mask_rna, col_sam_id]
        sams_dna_rna = list(set(sams_dna).intersect(set(sams_rna)))
        mask_qc = df_cln[col_sam_id].isin(sams_dna_rna)
    elif seqs_select.lower()=="all":
        mask_qc = mask_dna | mask_rna

    print("-selected %d/%d samples from seqs_select='%s'" % (sum(mask_qc), len(mask_qc), seqs_select))
    df_cln = df_cln.loc[mask_qc]

    # select subjects by cohort
    if subs_select!="all":
        if type(subs_select)==list:
            sub_list = subs_select
        else:
            if subs_select.lower()=="all":
                sub_list = df_cln[col_sub_id].tolist()
            elif subs_select.lower()=="breast":
                sub_list = df_cln.loc[df_cln["Cohort"].str.startswith("BREAST"), col_sub_id].tolist()
            elif subs_select.lower()=="lung":
                sub_list = df_cln.loc[df_cln["Cohort"].str.startswith("LUNG"), col_sub_id].tolist()

        mask_sub = df_cln[col_sub_id].isin(sub_list)
        df_cln = df_cln.loc[mask_sub].copy()
        print("-selected %d/%d samples from subs_select='%s'" % (sum(mask_sub), len(mask_sub), subs_select))

    # select samples
    dfs_cln = []
    subs_bas = df_cln.loc[df_cln["Cohort"].str.endswith("BAS"), col_sub_id].tolist()
    subs_ont = df_cln.loc[df_cln["Cohort"].str.endswith("ONT"), col_sub_id].tolist()
    subs_eot = df_cln.loc[df_cln["Cohort"].str.endswith("EOT"), col_sub_id].tolist()

    sams_bas = df_cln.loc[df_cln["Cohort"].str.endswith("BAS"), col_sam_id].tolist()
    sams_ont = df_cln.loc[df_cln["Cohort"].str.endswith("ONT"), col_sam_id].tolist()
    sams_eot = df_cln.loc[df_cln["Cohort"].str.endswith("EOT"), col_sam_id].tolist()

    if sams_select=="all":
        dfs_cln.append(df_cln)
    if "bas" in sams_select.split("_"):
        dfs_cln.append(df_cln.loc[df_cln[col_sam_id].isin(sams_bas)].copy())
    if "ont" in sams_select.split("_"):
        dfs_cln.append(df_cln.loc[df_cln[col_sam_id].isin(sams_ont)].copy())
    if "eot" in sams_select.split("_"):
        dfs_cln.append(df_cln.loc[df_cln[col_sam_id].isin(sams_eot)].copy())
    elif "baseot" in sams_select.split("_"):
        subs_baseot = list(set(subs_bas).intersection(set(subs_eot)))
        dfs_cln.append(df_cln.loc[df_cln[col_sub_id].isin(subs_baseot)].copy())
    elif "onteot" in sams_select.split("_"):
        subs_onteot = list(set(subs_ont).intersection(set(subs_eot)))
        dfs_cln.append(df_cln.loc[df_cln[col_sub_id].isin(subs_onteot)].copy())
    elif "basonteot":
        subs_basonteot = list(set(set(subs_bas).intersection(set(subs_eot))).intersection(set(subs_ont)))
        dfs_cln.append(df_cln.loc[df_cln[col_sub_id].isin(subs_basonteot)].copy())
    n_bef = df_cln.shape[0]
    df_cln = pd.concat(dfs_cln)
    n_aft = df_cln.shape[0]
    print("-selected %d/%d samples from sams_select='%s'" % (n_aft, n_bef, sams_select))

    # final sample list
    sam_list = df_cln["Sample_Id"].dropna().tolist()

    # process alterations
    df_alt = process_alt(df_alt, df_cna, df_mut, sam_list, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                         keep_alt_det, keep_alt_det_one, mut_det_level)

    return df_alt, df_cln


def group_colocalized_cnas(df_alt, df_hgnc, df_gff3, col_alt, col_gen, col_alt_cat, col_alt_det, cols_cln):
    df_hgnc_use = df_hgnc.loc[df_hgnc["Status"]=="Approved"].copy()
    df_hgnc_use = df_hgnc_use.loc[~df_hgnc["Chromosome"].isnull()].copy()

    # some alterations affect colocalized genes, treat them as one event
    alts_groups = ["1q32.1 Amp", "9p21.3 Del", "11q13 Amp", "11q14 Amp", "19p13 Amp", "19p13.3 Del", "19q13.32 Del",
                   "8p11.23 Amp", "20q13.33 Amp", "20q13.33 Del", "8p11.22 Amp" , "17q23 Amp", "22q11.2 Del"]
    print("-grouping colocalized CNA in the following locations:")
    print("\t" + "\n\t".join(alts_groups))
    cna_long2short={"Amp": ["Amplification", "HL Amp", "ML Amp", "LL Amp"],
                    "Del": ["Deletion", "Hom Del", "LOH"]}

    # modify df_alt if necessary to add "Del/Amp" suffix in alterations not having this suffix
    for alts_group in alts_groups:
        loc = alts_group.split()[0]
        cna = alts_group.split()[-1]
        genes_group = df_hgnc_use.loc[df_hgnc_use["Chromosome"].str.startswith(loc), "Approved symbol"].tolist()
        genes_group += df_gff3.loc[df_gff3["Chromosome"].str.startswith(loc), "Symbol"].tolist()
        mask_cna = df_alt[col_alt_cat].isin(cna_long2short[cna])
        genes_group = sorted(list(set(genes_group).intersection(set(df_alt.loc[mask_cna, col_gen]))))
        if len(genes_group)<=1:
            print("WARNING! %s has less than 2 members: %s" % (alts_group, genes_group))
        else:
            mask_genes = df_alt[col_gen].isin(genes_group)
            df_alt.loc[mask_genes & mask_cna, col_alt_det] = \
                    df_alt.loc[mask_genes & mask_cna, [col_gen, col_alt_det]].apply(" ".join, axis=1)
            df_alt.loc[mask_genes & mask_cna, col_gen] = alts_group
            df_alt.loc[mask_genes & mask_cna, col_alt] = alts_group

    # group rows if multiple genes affected in the same sample
    df_alt["Order"] = np.arange(df_alt.shape[0])
    cna_alt_cats = [e for x in cna_long2short.values() for e in x]
    mask_cna_gby = df_alt[col_alt_cat].isin(cna_alt_cats) & df_alt[col_alt].isin(alts_groups)
    df_alt_cna = df_alt.loc[mask_cna_gby].copy()
    df_alt_oth = df_alt.loc[~mask_cna_gby].copy()

    cols_gby = [x for x in df_alt_cna if x.endswith("_Id") if not "Civic" in x or "Oncokb" in x]
    cols_gby += [x for x in cols_cln if x in df_alt_cna and not x in cols_gby]
    cols_gby += [col_gen, col_alt]
    cols_oth = [x for x in df_alt_cna if x not in cols_gby+["Order", col_alt_cat]]
    agg_func = lambda  x: " & " .join(x.fillna("N/A").astype(str).tolist()) if not x.isnull().mean()==1 else np.nan
    df_alt_cna[cols_gby] = df_alt_cna[cols_gby].fillna("N/A")
    dt_agg_oth = {**{col_oth:  agg_func for col_oth in cols_oth}, **{"Order": "min"}, **{col_alt_cat: "first"}}
    df_alt_cna = df_alt_cna.groupby(cols_gby).agg(dt_agg_oth).reset_index()
    df_alt_cna[cols_gby] = df_alt_cna[cols_gby].replace({"N/A": np.nan})
    df_alt = pd.concat((df_alt_cna, df_alt_oth)).sort_values(by="Order")
    del df_alt["Order"]

    return df_alt


def group_by_pathway(df_alt, df_pth, col_alt, col_gen, col_alt_cat, col_alt_det, col_t_vaf, cols_cln):
    # group alterations by pathway
    alts_groups = df_pth["Pathway_Name"].drop_duplicates().tolist()

    # modify df_alt if necessary to add "Del/Amp" suffix in alterations not having this suffix
    for alts_group in alts_groups:
        genes_group = df_pth.loc[df_pth["Pathway_Name"]==alts_group, "Symbol"].tolist()
        genes_group = sorted(list(set(genes_group).intersection(set(df_alt[col_gen]))))
        if len(genes_group)<=1:
            print("WARNING! %s has less than 2 members: %s" % (alts_group, genes_group))

        mask_genes = df_alt[col_gen].isin(genes_group)
        df_alt.loc[mask_genes, col_alt_det] = \
                df_alt.loc[mask_genes, [col_gen, col_alt_det]].apply(" ".join, axis=1)
        df_alt.loc[mask_genes, col_gen] = alts_group
        df_alt.loc[mask_genes, col_alt] = alts_group

    # group rows if multiple genes affected in the same sample
    df_alt["Order"] = np.arange(df_alt.shape[0])

    # keep detail of col_alt_cat per aggregated alteration
    df_alt[f"{col_alt_cat}_Detail"] = df_alt[col_alt_cat]

    cols_gby = [x for x in df_alt if x.endswith("_Id") if not "Civic" in x or "Oncokb" in x]
    cols_gby += [x for x in cols_cln if x in df_alt and not x in cols_gby]
    cols_gby += [col_gen, col_alt]
    cols_con = [col_alt_cat, "Clinical_Benefit", "Responder"]
    cols_con = [x for x in cols_con if x in df_alt and x not in cols_gby]
    cols_oth = [x for x in df_alt if x not in cols_gby+cols_con+["Order", col_t_vaf]]
    agg_func_gen = lambda  x: " & " .join(x.fillna("N/A").astype(str).tolist()) if not x.isnull().mean()==1 else np.nan
    agg_func_con = lambda  x: " & " .join(sorted(x.drop_duplicates().fillna("N/A").astype(str).tolist()))
    agg_func_vaf = lambda  x: max(x) if x.isnull().sum()==0 else np.nan

    df_alt[cols_gby] = df_alt[cols_gby].fillna("N/A")
    dt_agg_oth = {**{col_oth: agg_func_gen for col_oth in cols_oth}, **{"Order": "min"},
                  **{col_con: agg_func_con for col_con in cols_con}, **{col_t_vaf: agg_func_vaf}}
    df_alt = df_alt.groupby(cols_gby).agg(dt_agg_oth).reset_index()
    df_alt[cols_gby] = df_alt[cols_gby].replace({"N/A": np.nan})
    df_alt = df_alt.sort_values(by="Order")
    del df_alt["Order"]

    return df_alt

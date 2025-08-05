# -*- coding: utf-8 -*-
"""
@created: May 03 2022
@modified: Sep 09 2024
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Oncoplot-like figure detailing treatment resistances
"""

import argparse
import numpy as np
import pandas as pd
import sys
sys.path.append("../common/functions")

from utils import combine_all_alterations
from utils_comut import add_responder_status

# functions ============================================================================================================

def main(args):
    # define util variable names
    col_gen = "Hugo_Symbol"
    col_alt = "Alteration"
    col_alt_det = "Alteration_Detail"
    col_sub_id = "Subject_Id"
    col_sam_id = "Sample_Id"
    col_alt_cat = "Alteration_Category"

    df_alt, df_cln = combine_all_alterations(alt=args.alt, cna=args.cna, mut=args.mut, cln=args.cln,
                                             col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                             col_sub_id=col_sub_id, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                             keep_alt_det=True, cna_selection="focal_non_neutral",
                                             non_syn_mut_only=False, seqs_select="dna", subs_select="all",
                                             sams_select="all")

    # select data for gen and add rows for samples not altered
    cols_cln = [col_sub_id, col_sam_id]
    cols_alt = [col_sub_id, col_sam_id, col_gen, col_alt, col_alt_det, "t_vaf", "Cellular_Fraction", "Annotated"]
    df_alt_gen = df_alt.loc[df_alt[col_gen]==args.gen]
    df_fin_gen = df_cln[cols_cln].merge(df_alt_gen[cols_alt], how="left", on=[col_sub_id, col_sam_id])
    df_fin_gen = df_fin_gen.sort_values(by=[col_sub_id, col_sam_id])
    mask_nna = ~df_fin_gen["t_vaf"].isnull()
    df_fin_gen.loc[mask_nna, "t_vaf"] = df_fin_gen.loc[mask_nna, "t_vaf"].apply(lambda x: "%.3g %%" % (x*100))
    mask_nna = ~df_fin_gen["Cellular_Fraction"].isnull()
    df_fin_gen.loc[mask_nna, "Cellular_Fraction"] = \
            df_fin_gen.loc[mask_nna, "Cellular_Fraction"].apply(lambda x: "%.3g %%" % (x*100))

    # empty col_alt when equals to gene symbol
    mask_gen = df_fin_gen[col_alt]==df_fin_gen[col_gen]
    df_fin_gen.loc[mask_gen, col_alt] = np.nan

    # add clinical columns
    cols_cln = ["Subject_Id", "OS_Time_Days", "OS_Status",  "PFS_Time_Days", "PFS_Status", "Best_Overall_Response",
                "Clinical_Benefit", "Cancer_Type"]
    df_cln = add_responder_status(df_cln)
    df_fin_gen = df_fin_gen.merge(df_cln[cols_cln].drop_duplicates(), how="left", on="Subject_Id")
    df_fin_gen["Timepoint"] = df_fin_gen["Sample_Id"].apply(lambda x: x.split("-")[-2])

    # create table with wild-type or altered status
    cols_gby = ["Sample_Id"]
    cols_agg = [x for x in df_fin_gen if x not in cols_gby]
    dt_agg = {**{"Alteration": lambda x: ";".join(x.dropna())}, **{x: "first" for x in cols_agg if x!="Alteration"}}
    df_gen = df_fin_gen.groupby(cols_gby).agg(dt_agg)
    df_gen["Status"] = "Wild-type"
    maks_alt = df_gen["Alteration"]!=""
    df_gen.loc[maks_alt, "Status"] = "Altered"

    ## cross tables
    df_gen_bas = df_gen.loc[df_gen["Timepoint"]=="BAS"]

    # Clinical Benefit
    df_cnt_cb = df_gen_bas.groupby(["Cancer_Type", "Status", "Clinical_Benefit"]).size()
    df_cnt_cb = df_cnt_cb.unstack(level=[0,1])
    df_cnt_cb.index = pd.MultiIndex.from_product([[""], df_cnt_cb.index], names=["",""])

    # BOR
    df_cnt_bor = df_gen_bas.groupby(["Cancer_Type", "Status", "Best_Overall_Response"]).size()
    df_cnt_bor = df_cnt_bor.unstack(level=[0,1])
    df_cnt_bor.index = pd.MultiIndex.from_product([[""], df_cnt_bor.index], names=["", ""])

    # Timepoint
    df_cnt_tp = df_gen.groupby(["Cancer_Type", "Status", "Timepoint"]).size()
    df_cnt_tp = df_cnt_tp.unstack(level=[0,1])
    df_cnt_tp.index = pd.MultiIndex.from_product([[""], df_cnt_tp.index], names=["",""])

    # aggregate into one
    multi_index_cb  = pd.MultiIndex.from_tuples([("Clinical_Benefit (BAS samples)", "")], names=["", ""])
    multi_index_bor  = pd.MultiIndex.from_tuples([("Best Objective Response (BAS samples)", "")], names=["", ""])
    multi_index_tp  = pd.MultiIndex.from_tuples([("Timepoint", "")], names=["", ""])


    df_cnt_cb = pd.concat((pd.DataFrame(index=multi_index_cb), df_cnt_cb, pd.DataFrame(index=multi_index_bor)))
    df_cnt_bor = pd.concat((df_cnt_bor, pd.DataFrame(index=multi_index_tp)))
    df_cnt = pd.concat((df_cnt_cb, df_cnt_bor, df_cnt_tp), axis=0)

    # select patients
    if args.select:
        df_sub = df_fin_gen[[col_sub_id, "t_vaf"]].dropna(how="any")
        df_fin_gen = df_fin_gen.loc[df_fin_gen[col_sub_id].isin(df_sub[col_sub_id])]

    # save
    with pd.ExcelWriter(args.output) as writer:
        df_fin_gen.to_excel(writer, index=False, sheet_name="Alterations", float_format="%.2f")
        df_cnt.to_excel(writer, index=True, sheet_name="Summaries")

    print("-file saved at %s" % args.output)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Oncoplot-like figure detailing treatment resistances.")
    parser.add_argument("--cln", type=str, help="Path to input curated clincal table.",
                        default="../../data/cln/curated/cln_icarus_curated_with_qc.tsv")
    parser.add_argument('--alt', type=str, help='Path to table of alterations.',
                        default="../../results/wes_analysis/alterations/aggregated_alterations.tsv")
    parser.add_argument('--mut', type=str, help='Path to table of all mutations.',
                        default="../../data/wes/somatic_maf/somatic_calls.maf.gz")
    parser.add_argument('--cna', type=str, help='Path to table of all CNAs.',
                        default="../../data/wes/somatic_cna/somatic_calls.tsv.gz")
    parser.add_argument('--gen', type=str, help='Gene symbol', default="TACSTD2")
    parser.add_argument('--select', action="store_true", default=False,
                        help='If True, select only samples from patients bearing an alteration in at least 1 sample')
    parser.add_argument('--output', type=str,  help='Paths to output table.',
                        default="../../results/wes_analysis/alterations/alterations_TACSTD2_all.xlsx")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)

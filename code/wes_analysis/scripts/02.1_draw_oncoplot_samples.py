# -*- coding: utf-8 -*-
"""
@created: May 02 2022
@modified: Aug 05 2025
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Oncoplot-like figure.
"""

import argparse
import numpy as np
import os
import pandas as pd
import re
import sys

import matplotlib.pyplot as plt

# utils
sys.path.append("../common/functions")
from utils import combine_all_alterations, group_colocalized_cnas, group_by_pathway
from utils_comut import get_ordered_alterations_and_samples, add_vaf_for_mutations, add_responder_status
from utils_comut import get_tables, get_mappings, draw_comut_plot, run_fisher_boschloo_tests
from utils_comut import select_alt_from_gene_lists, select_inactivating_alt

def anonymize_ids(sam, ids_anonym, col_sub_id, col_sub_anonym_id, col_sam_id, col_sam_anonym_id, df_alt, df_cln):
    # use Sample Id and anonymize
    df_sam = pd.read_table(sam)
    df_anonym = pd.read_table(ids_anonym)
    df_sam = df_sam.merge(df_anonym, how="left")

    # create sample suffix from sample id. Serve to remove biopsy number in rare cases of double biopsies at the same
    # timepoint
    col_bio_id = "Biopsy_Id"
    def match_pattern(x, pattern):
        match = re.search(pattern, x)
        if match:
            return match.group(0)
        else:
            return ""
    regex_sam_suffix = r"-?(C[1-9]J[0-9]{1,2}|BAS|EOT)?-?[TN]"
    df_sam["Sample_Suffix"] = df_sam[col_bio_id].apply(match_pattern, pattern=regex_sam_suffix)
    df_sam[col_sam_anonym_id] = df_sam[col_sub_anonym_id] + df_sam["Sample_Suffix"]

    anonym_sub_ids = {r[col_sub_id]: r[col_sub_anonym_id] for _, r in \
        df_sam[[col_sub_id, col_sub_anonym_id]].dropna(subset=[col_sub_id]).drop_duplicates().iterrows()}
    anonym_sam_ids = {r[col_bio_id]: r[col_sam_anonym_id] for _, r in \
        df_sam[[col_bio_id, col_sam_anonym_id]].dropna(subset=[col_bio_id]).drop_duplicates().iterrows()}

    # anonymise ids
    for df in [df_alt, df_cln]:
        df[col_sub_id] = df[col_sub_id].replace(anonym_sub_ids)
        df[col_sam_id] = df[col_sam_id].replace(anonym_sam_ids)

    return df_alt, df_cln


def find_overlapping_bands(gene_chr, gene_start, gene_end, bands):
    bands_on_same_chr = bands[bands['Chromosome'] == gene_chr]
    overlapping_bands = bands_on_same_chr[(bands_on_same_chr['Start'] <= gene_end) & \
                                          (bands_on_same_chr['End'] >= gene_start)]
    return ';'.join(overlapping_bands['Band'].tolist())


def load_gff3_and_add_bands(gff3, band):
    df_gff3 = pd.read_table(gff3, header=None)
    df_gff3.columns = ["Chromosome", "Start", "End", "Ensembl_Id", "Symbol", "Biotype", "Source"]
    df_band = pd.read_table(band, header=None)
    df_band.columns = ["Chromosome", "Start", "End", "Band", "Other"]
    df_band["Chromosome"] = df_band["Chromosome"].apply(lambda x: re.sub("^chr", "", x))
    df_gff3["Band"] = \
        df_gff3.apply(lambda row: find_overlapping_bands(row['Chromosome'], row['Start'], row['End'], df_band), axis=1)
    df_gff3["Band"] = df_gff3["Band"].str.split(";")
    df_gff3 = df_gff3.explode("Band")
    df_gff3["Chromosome"] = df_gff3["Chromosome"] + df_gff3["Band"]
    return df_gff3


def main(args):
    # define util variable names
    col_gen = "Hugo_Symbol"
    col_alt = "Alteration"
    col_alt_cat = "Alteration_Category"
    col_alt_cat_sec = "Alteration_Category_Secondary"
    col_alt_det = "Alteration_Detail"
    col_sub_id = "Subject_Id"
    col_sam_id = "Sample_Id"
    col_t_vaf = "t_vaf"
    font_mode = "not_tex"
    col_sub_anonym_id = "%s_Anonymous" % col_sub_id
    col_sam_anonym_id = "%s_Anonymous" % col_sam_id

    if font_mode=="tex":
        plt.rcParams['text.usetex'] = True
    else:
        plt.rcParams['text.usetex'] = False

    # load data
    if args.annotated_only=="yes":
        df_alt, df_cln = combine_all_alterations(alt=args.alt, cln=args.cln, cna=None, mut=None, col_gen=col_gen,
                                                 col_alt=col_alt, col_alt_det=col_alt_det, col_alt_cat=col_alt_cat,
                                                 col_sub_id=col_sub_id, col_sam_id=col_sam_id,
                                                 seqs_select=args.seqs_select, subs_select=args.subs_select,
                                                 sams_select=args.sams_select, keep_alt_det=False,
                                                 cna_selection="focal_high_level")

        if args.extended=="yes":
            df_alt_ext, _ = combine_all_alterations(alt=args.alt, cln=args.cln, cna=args.cna, mut=args.mut,
                                                    col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                                    col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                                                    col_sam_id=col_sam_id, seqs_select=args.seqs_select,
                                                    subs_select=args.subs_select, sams_select=args.sams_select,
                                                    keep_alt_det=False, keep_alt_det_one=False, non_syn_mut_only=True,
                                                    cna_selection="focal_non_neutral")

            genes_ann = df_alt[col_gen].unique().tolist()
            mask_gen = df_alt_ext[col_gen].isin(genes_ann)
            df_alt = df_alt_ext.loc[mask_gen].copy()
            print("-INFO: selected %d/%d extended alterations from a list of %d genes with drivers" % \
                  (sum(mask_gen), len(mask_gen), len(genes_ann)))


    else:
        if args.extended=="yes":
            df_alt, df_cln = combine_all_alterations(alt=args.alt, cln=args.cln, cna=args.cna, mut=args.mut,
                                                     col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                                     col_alt_cat=col_alt_cat, col_sub_id=col_sub_id, col_sam_id=col_sam_id,
                                                     seqs_select=args.seqs_select, subs_select=args.subs_select,
                                                     sams_select=args.sams_select, keep_alt_det=False,
                                                     keep_alt_det_one=False, non_syn_mut_only=True,
                                                     cna_selection="focal_non_neutral")
        else:
            df_alt, df_cln = combine_all_alterations(alt=args.alt, cln=args.cln, cna=args.cna, mut=args.mut,
                                                     col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                                     col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                                                     col_sam_id=col_sam_id, seqs_select=args.seqs_select,
                                                     subs_select=args.subs_select, sams_select=args.sams_select,
                                                     keep_alt_det=False, non_syn_mut_only=True,
                                                     cna_selection="focal_high_level")


    # rename/reformat some categorical values
    df_cln["Cancer_Type"] = df_cln["Cancer_Type"].replace({"Breast Cancer": "BREAST", "Lung Cancer": "LUNG"})
    df_cln["BOR_Confirmed"] = df_cln["BOR_Confirmed"].fillna("N/A")
    df_cln["Best_Overall_Response"] = df_cln["Best_Overall_Response"].fillna("N/A")
    df_cln["MSKCC_Oncotree"] = df_cln["MSKCC_Oncotree"].replace({"Lung": "N/A", "BRCA": "N/A"}).fillna("N/A")
    df_cln["Clinical_Benefit"] = df_cln["Clinical_Benefit"].fillna("N/A")

    # add binary benefit for plot
    df_cln = add_responder_status(df_cln)

    # add Match_DNA column
    mask_matched = ~df_cln["Sample_Id_DNA_N"].isnull()
    df_cln.loc[mask_matched, "Matched_DNA"] = "Matched normal"
    df_cln.loc[~mask_matched, "Matched_DNA"] = "Unmatched"

    # remap Cancer_Type
    df_cln["MSKCC_Oncotree"] = df_cln["MSKCC_Oncotree"].replace({"LUAD": "NSQ", "LUAS": "SCC", "LUSC": "SCC"})
    mask_nsq = df_cln["Diagnosis_Histology"]=="Non-squamous"
    df_cln.loc[mask_nsq, "MSKCC_Oncotree"] = "NSQ"

    # add/harmonize clinical columns
    cols_cln_req = ["Cohort", "Timepoint", "Cancer_Type", "Clinical_Benefit", "Responder",
                    "BOR_Confirmed", "Best_Overall_Response", "MSKCC_Oncotree"]
    cols_cln_req += [x for x in df_alt if x in df_cln and x not in cols_cln_req + [col_sub_id, col_sam_id]]

    for col_cln in cols_cln_req:
        if col_cln in df_alt:
            del df_alt[col_cln]
        df_alt = df_alt.merge(df_cln[[col_sam_id, col_cln]], how="left", on=col_sam_id)
        assert df_alt[col_cln].isnull().sum()==0

    cols_cln = df_cln.columns.tolist()

    # keep intact alterations table for top barplot
    if args.top_mode=="all":
        df_alt_top = df_alt.copy()
    else:
        df_alt_top = None

    # restrict to alterations from specified gene lists
    if args.gene_list_names is not None and len(args.gene_list_names)>0:
        df_alt = select_alt_from_gene_lists(df_alt=df_alt, col_gen=col_gen, gene_list_path=args.gene_list_path,
                                            gene_list_names=args.gene_list_names)

    # restrict alteration description to gene name
    if args.gene_only=="yes":
        df_alt[col_alt] = df_alt[col_gen]


    # restrict to inactivating alterations
    if args.inactivating_only=="yes":
        df_alt = select_inactivating_alt(df_alt=df_alt, col_alt_cat=col_alt_cat)

    # group by pathway
    df_hgnc = pd.read_table(args.hgnc)
    df_gff3 = load_gff3_and_add_bands(args.gff3, args.band)

    if args.pathway_use=="yes":
        n_alt_bef = df_alt.shape[0]

        # add pathway tag to avoid mixing pathway names with gene names
        df_pth = pd.read_table(args.pathway_path)
        df_pth["Pathway_Name"] = df_pth["Pathway_Name"] + " Pathway"
        df_alt = group_by_pathway(df_alt, df_pth, col_alt, col_gen, col_alt_cat, col_alt_det, col_t_vaf, cols_cln)

        # retain only pathways
        mask_pth = df_alt[col_alt].str.endswith(" Pathway")
        df_alt = df_alt.loc[mask_pth].copy()
        print("-INFO: dropping %d/%d alterations not in any of the pathways provided" % (sum(~mask_pth), n_alt_bef))

        # remove suffix
        df_alt[col_alt] = df_alt[col_alt].apply(lambda x: re.sub(" Pathway$", "", x))
        df_pth["Pathway_Name"] = df_pth["Pathway_Name"].apply(lambda x: re.sub(" Pathway$", "", x))
    else:
        # some alterations affect colocalized genes, treat them as one event
        df_alt =  group_colocalized_cnas(df_alt, df_hgnc, df_gff3, col_alt, col_gen, col_alt_cat, col_alt_det, cols_cln)

    # data for top bar plot
    if args.top_mode=="all":
        # restrict to inactivating alterations
        if args.inactivating_only=="yes":
            df_alt_top = select_inactivating_alt(df_alt=df_alt_top, col_alt_cat=col_alt_cat)

        # grouped alterations affecting colocalized genes in the counting of events on top bar plot
        df_alt_top =  group_colocalized_cnas(df_alt_top, df_hgnc, df_gff3, col_alt, col_gen, col_alt_cat, col_alt_det,
                                             cols_cln)
    else:
        if args.pathway_use=="no":
            df_alt_top = df_alt.copy()
        else:
            df_alt_top = df_alt[[col_sam_id, col_alt, col_alt_det, f"{col_alt_cat}_Detail"]]
            df_alt_top = df_alt_top.rename(columns={f"{col_alt_cat}_Detail": col_alt_cat})

            df_alt_top[col_alt_det] = df_alt_top[col_alt_det].str.split(' & ')
            df_alt_top[col_alt_cat] = df_alt_top[col_alt_cat].str.split(' & ')
            df_alt_top = df_alt_top.explode([col_alt_det, col_alt_cat], ignore_index=True)


    # if alterations were grouped by pathways, make a secondary alteration column
    if args.pathway_use=="yes":
        # simplify categories to retain only the first and secondary, if any
        extract_alt_cat_sec = lambda x: x.split(" & ")[1] if " & " in x else np.nan
        remove_alt_cat_sec = lambda x: x.split(" & ")[0]
        df_alt[col_alt_cat_sec] = df_alt[col_alt_cat].apply(extract_alt_cat_sec)
        df_alt[col_alt_cat] = df_alt[col_alt_cat].apply(remove_alt_cat_sec)

        # if at least one alteration in the pathway is annotated, then show indicator
        mask_ann = df_alt["Annotated"].str.contains("Yes")
        df_alt.loc[mask_ann, "Annotated"] = "Yes"
        df_alt.loc[~mask_ann, "Annotated"] = "No"
    else:
        df_alt[col_alt_cat_sec] = np.nan

    # simplify categories
    alt_cat_simple = {"Amplification": "HL Amp", "cnLOH": "LOH", "LL_gain": "LL Amp", "ML_gain": "ML Amp",
                      "Deletion": "Hom Del"}
    df_alt[col_alt_cat] = df_alt[col_alt_cat].replace(alt_cat_simple)
    df_alt_top[col_alt_cat] = df_alt_top[col_alt_cat].replace(alt_cat_simple)
    df_alt[col_alt_cat_sec] = df_alt[col_alt_cat_sec].replace(alt_cat_simple)

    # use Sample Id and anonymize
    if args.anonymize.lower()=="yes":
        df_alt, df_cln = anonymize_ids(args.sam, args.ids_anonym, col_sub_id, col_sub_anonym_id, col_sam_id,
                                       col_sam_anonym_id, df_alt, df_cln)
        df_alt_top, _ = anonymize_ids(args.sam, args.ids_anonym, col_sub_id, col_sub_anonym_id, col_sam_id,
                                      col_sam_anonym_id, df_alt_top, df_cln)

    # parameters for the plot
    sams = df_cln[col_sam_id].tolist()
    cols_comut = ["sample", "category", "value"]

    if args.sams_select=="baseot":
        threshold_alt_sams = df_cln.loc[df_cln["Timepoint"]=="EOT", col_sam_id].tolist()
    else:
        threshold_alt_sams = None

    alts_ordered, sams_ordered = get_ordered_alterations_and_samples(df_cln=df_cln, df_alt=df_alt, sams=sams,
                                                                     col_gen=col_gen, col_alt=col_alt,
                                                                     col_sub_id=col_sub_id, col_sam_id=col_sam_id,
                                                                     threshold_alt=args.threshold_alt,
                                                                     threshold_alt_sams=threshold_alt_sams,
                                                                     complete_gene=True, contiguous_sub=True)


    # identify the variable underlying the left barplot
    if args.col_side.lower()=="benefit":
        col_side_barplot = "Clinical_Benefit"
        name_side_barplot = "Clinical Benefit"
        values_side_barplot_ignore = ["N/A", "NE"]
    elif args.col_side.lower()=="response":
        col_side_barplot = "Responder"
        name_side_barplot = "Responder"
        values_side_barplot_ignore = ["N/A", "NE"]
    elif args.col_side.lower()=="timepoint":
        col_side_barplot = "Timepoint"
        name_side_barplot = "Timepoint"
        values_side_barplot_ignore = []
    else:
        raise ValueError(f"Value '{args.col_side}' for --col_side is not supported. Choose 'benefit', 'response'" + \
                         " or 'timepoint'")


    # order samples based on clinical features
    if col_side_barplot in ["Clinical_Benefit", "Responder"]:
        if col_side_barplot=="Clinical_Benefit":
            cols_sort = ["Clinical_Benefit", "Responder", "Best_Overall_Response", "Order"]
        elif col_side_barplot=="Responder":
            cols_sort = ["Responder", "Best_Overall_Response", "Order"]
    else:
        cols_sort = ["Best_Overall_Response", "Order"]

    df_cln["Responder"] = pd.Categorical(df_cln["Responder"], ["Yes", "No"])
    df_cln["Clinical_Benefit"] = pd.Categorical(df_cln["Clinical_Benefit"], ["Yes", "No"])
    categories_response = ["CR", "PR", "PR conf.", "PR unconf.", "SD", "PD", "NE", "N/A"]
    categories_response = [x for x in categories_response if x in df_cln["Best_Overall_Response"].unique()]
    df_cln["Best_Overall_Response"] = pd.Categorical(df_cln["Best_Overall_Response"], categories_response)
    df_alt["Best_Overall_Response"] = pd.Categorical(df_alt["Best_Overall_Response"], categories_response)
    df_sam_order = pd.DataFrame({col_sam_id: sams_ordered, "Order": np.arange(0, len(sams_ordered))})
    df_cln = df_cln.merge(df_sam_order, how="left", on=col_sam_id)

    if args.subs_select=="all":
        cols_sort = ["Cancer_Type"] + cols_sort

    df_cln = df_cln.sort_values(by=cols_sort)
    sams_ordered = df_cln[col_sam_id].tolist()

    # if not using grouping by pathways,
    #   add rows with VAF category for col_alt_cat in order to represent the VAF in the oncoplot
    # else
    #   add rows with col_alt_cat_sec if any, or VAF category otherwise
    df_alt, t_vaf_labs = add_vaf_for_mutations(df_alt=df_alt, df_cln=df_cln, col_sub_id=col_sub_id,
                                               col_sam_id=col_sam_id, col_alt=col_alt, col_alt_det=col_alt_det,
                                               col_alt_cat=col_alt_cat, col_t_vaf=col_t_vaf,
                                               col_alt_cat_sec=col_alt_cat_sec, mode=font_mode)


    # data for the plot

    if args.subs_select=="all":
        categorical_data = ["tim", "his", "can", "res", "mtc"]
        categorical_names = ["Timepoint", "Histology",  "Cancer type", "Response", "Matched DNA"]
        categorical_vars = ["Timepoint", "MSKCC_Oncotree", "Cancer_Type", "Best_Overall_Response",
                            "Matched_DNA"]
    else:
        if col_side_barplot=="Responder":
            categorical_data = ["res", "res_bin", "mtc"]
            categorical_names = ["Best Response", "Responder", "Matched DNA"]
            categorical_vars = ["Best_Overall_Response", "Responder", "Matched_DNA"]
        elif col_side_barplot=="Clinical_Benefit":
            categorical_data = ["res", "res_bin", "ben", "mtc"]
            categorical_names = ["Best Response", "Responder", "Clinical Benefit", "Matched DNA"]
            categorical_vars = ["Best_Overall_Response", "Responder", "Clinical_Benefit", "Matched_DNA"]
        elif col_side_barplot=="Timepoint":
            categorical_data = ["tim", "res", "res_bin", "ben", "mtc"]
            categorical_names = ["Timepoint", "Best Response", "Responder", "Clinical Benefit", "Matched DNA"]
            categorical_vars = ["Timepoint", "Best_Overall_Response", "Responder", "Clinical_Benefit", "Matched_DNA"]

        if args.subs_select=="lung":
            categorical_data.insert(1, "his")
            categorical_names.insert(1, "Histology")
            categorical_vars.insert(1, "MSKCC_Oncotree")

    categorical_dnv = {d: {"name": n, "var": v} for d,n,v in zip(categorical_data, categorical_names, categorical_vars)}

    # get tables and mappings for the plot
    dfs_data = get_tables(df_cln=df_cln, df_alt=df_alt, col_alt=col_alt, col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                          col_sam_id=col_sam_id, alts_ordered=alts_ordered, categorical_dnv=categorical_dnv,
                          col_side_barplot=col_side_barplot, values_side_barplot_ignore=values_side_barplot_ignore,
                          df_alt_top=df_alt_top)
    mappings = get_mappings(t_vaf_labs=t_vaf_labs, col_side_barplot=col_side_barplot)

    height_middle = 0.1*len(alts_ordered)
    ncol_legend = 1
    if len(alts_ordered)<=15:
        height_top = 6
        if len(alts_ordered) < 8:
            bbox_to_anchor_legend = (1,2)
        else:
            bbox_to_anchor_legend = (1,1.3)
        shadow_width_left = 1.5
        ytick_fontsize = 12
        ncol_legend = 2
    elif len(alts_ordered)<=30:
        height_top = 6
        bbox_to_anchor_legend = (1,1)
        shadow_width_left = 1.75
        ytick_fontsize = 12
        ncol_legend = 2
    else:
        height_top = 8
        if len(sams_ordered) > 30:
            shadow_width_left = 1.2
        elif len(sams_ordered) > 20:
            shadow_width_left = 1.4
        else:
            shadow_width_left = 1.6
        bbox_to_anchor_legend = (1,1)
        ytick_fontsize = 10

    # if only gene names, can reduce
    if args.gene_only=="yes":
        shadow_width_left -= 0.3

    ignored_plots = ["Alteration Type", "Same patient bis"]
    ignored_values = ["Absent"] + t_vaf_labs

    labels_orders = {"Alterations": ["Indel", "Mutation", "HL Amp", "ML Amp", "LL Amp", "Hom Del", "LOH"] +
                                    ["Multiple", "Absent"] + t_vaf_labs}
    for d in ["res", "res_bin", "ben"]:
        if d in categorical_dnv:
            if d=="res":
                labels_orders[categorical_dnv[d]["name"]] = categories_response
            else:
                labels_orders[categorical_dnv[d]["name"]] = ["Yes", "No"]

    # draw figure
    if args.sams_select in ["all", "baseot"]:
        sample_indicators = True
        stacked_bar_side = False
        title_legend_barplot = col_side_barplot.replace("_", " ")
    else:
        sample_indicators = False
        stacked_bar_side = False
        title_legend_barplot = col_side_barplot.replace("_", " ")

    if col_side_barplot in categorical_vars:
        show_legend_barplot = False
    else:
        show_legend_barplot = True

    if args.annotated_only=="yes":
        show_symbols = False
        label_bar_top = "Drivers"
        digits_bar_size_fontsize = 7
    else:
        show_symbols = True
        digits_bar_size_fontsize = 7
        if args.inactivating_only == "yes":
            label_bar_top = "Inactivating\nAlterations"
        else:
            label_bar_top = "Nonsyn.\nAlterations"

    comut_alt = draw_comut_plot(dfs_data=dfs_data, mappings=mappings, sams_ordered=sams_ordered,
                                alts_ordered=alts_ordered, t_vaf_labs=t_vaf_labs, categorical_data=categorical_data,
                                categorical_names=categorical_names, label_bar_top=label_bar_top,
                                shadow_width_left=shadow_width_left,
                                xtick_fontsize=8, ytick_fontsize=ytick_fontsize,
                                mode_bar_side="Percent", label_bar_side="Fraction of samples",
                                sample_indicators=sample_indicators, height_middle=height_middle, height_top=height_top,
                                stacked_bar_side=stacked_bar_side, labels_orders=labels_orders,
                                ignored_values=ignored_values, label_bar_top_fontsize=14, label_bar_side_fontsize=16,
                                digits_bar_size_fontsize=digits_bar_size_fontsize, ncol_legend=ncol_legend,
                                show_legend_barplot=show_legend_barplot, title_legend_barplot=title_legend_barplot,
                                unit_width=0.2, unit_height=0.2, bbox_to_anchor_legend=bbox_to_anchor_legend,
                                show_symbols=show_symbols, width_left=2.75)
    # save figure
    output_pdf_dir = os.path.dirname(args.output_pdf)
    os.makedirs(output_pdf_dir, exist_ok=True)
    fig = comut_alt.figure
    fig.savefig(args.output_pdf, dpi=300, bbox_inches='tight')
    print("-plot saved at %s" % args.output_pdf)

    # table of alterations underlying plot
    df_alt_inplot = df_alt.loc[df_alt[col_alt].isin(alts_ordered)]
    cols_tab = [col_sub_id, col_sam_id, col_gen, col_alt_cat, col_alt, col_alt_det, col_t_vaf, col_side_barplot] + \
            [x for x in df_alt_inplot if x.startswith("Oncokb")] + \
            [x for x in df_alt_inplot if x.startswith("Civic")]
    df_alt_tab = df_alt_inplot[cols_tab].copy()

    # drop artificially repeated rows to allow for representing multiple categories in one cell
    df_alt_tab = df_alt_tab.drop_duplicates(subset=[col_sub_id, col_sam_id, col_gen, col_alt, col_alt_det],
                                            keep="first")

    # if asked, run fisher-boschloo tests
    if args.run_tests=="yes" and col_side_barplot in ["Responder", "Clinical_Benefit", "Timepoint"]:
        df_alt_tests = run_fisher_boschloo_tests(df_cln=df_cln, df_alt=df_alt_tab,
                                                 col_side_barplot=col_side_barplot, tests_alt=args.tests_alt,
                                                 sams=sams, col_gen=col_gen, col_alt=col_alt, col_sub_id=col_sub_id,
                                                 col_sam_id=col_sam_id, threshold_alt_tests=args.threshold_alt,
                                                 alts_ordered=alts_ordered)
    else:
        df_alt_tests = None

    # save tables of alterations and, if not None, of tests
    if args.output_tab is not None:
        with pd.ExcelWriter(args.output_tab) as writer:
            df_alt_tab.to_excel(writer, sheet_name="Alterations", index=False)
            if df_alt_tests is not None:
                df_alt_tests.to_excel(writer, sheet_name="Tests_vs_%s" % col_side_barplot, index=False)


# run ==================================================================================================================

if __name__ == "__main__":
    R_FOLDER = "../../results/wes_analysis/oncoplot"

    parser = argparse.ArgumentParser(description="Oncoplot figure.")
    parser.add_argument("--cln", type=str, help="Path to input curated clincal table.",
                        default="../../data/cln/curated/cln_icarus_curated_with_qc.tsv")
    parser.add_argument('--alt', type=str, help='Path to table of alterations.',
                        default="../../results/wes_analysis/alterations/aggregated_alterations.tsv")
    parser.add_argument('--mut', type=str, help='Path to table of all mutations.',
                        default="../../data/wes/somatic_maf/somatic_calls.maf.gz")
    parser.add_argument('--cna', type=str, help='Path to table of all CNAs.',
                        default="../../data/wes/somatic_cna/somatic_calls.tsv.gz")
    parser.add_argument("--sam", type=str, help="Paths to table of all samples.",
                        default="../../data/cln/curated/sam_icarus_curated.tsv")
    parser.add_argument("--ids_anonym", type=str, help="Paths to input table of anonymised patient ids.",
                        default="../../results/tables_paper/table_anonymisation.tsv")
    parser.add_argument("--anonymize", type=str, help="Set to 'yes' to anonymize ids.",
                        default="no")
    parser.add_argument("--hgnc", type=str, help="Path to HGNC table linking gene names to cytobands.",
                        default="../../data/resources/hgnc_all_symbols_03012022.tsv")
    parser.add_argument("--gff3", type=str, help="Path to gff3 table of genes used for CNV calling.",
                        default="../../data/resources/Homo_sapiens.GRCh37.87.gff3.gene.bed")
    parser.add_argument("--band", type=str, help="Path to table of chromosomal band regions.",
                        default="../../data/resources/cytoBand_hg19_ucsc.txt")
    parser.add_argument('--annotated_only', type=str, default="yes",
                        help='Only alterations annotated by oncokb/civic are shown.')
    parser.add_argument('--inactivating_only', type=str, default="no",
                        help='Only inactivating events are shown.')
    parser.add_argument('--threshold_alt', type=float, default=5,
                        help='Only alterations with a frequency/count of at least this value will be displayed')
    parser.add_argument('--extended', type=str, default="no",
                        help='A custom selection of extended alterations among selected genes will be shown.')
    parser.add_argument('--gene_list_path', type=str,
                        default="../../data/resources/Gene_Lists_Of_Interest.xlsx",
                        help='Path to the table containing the gene lists to restrict the analysis to.')
    parser.add_argument('--gene_list_names', nargs="+", type=str,
                        # default=["LUNG-specific", "All"],
                        # default=["BREAST-specific", "All"],
                        # default=["Inactiv_Leads_To_Resistance"],
                        # default=["Inactiv_Leads_To_Sensitivity"],
                        default=[],
                        help='Names of the gene lists to restrict the analysis to.')
    parser.add_argument("--pathway_path", type=str, help="Path to table of pathways.",
                        default="../../data/resources/pathways/data/sanchez-vega_cell_2018_curated.tsv")
    parser.add_argument("--pathway_use", type=str, default="no")
    parser.add_argument('--seqs_select', type=str, default="dna", help='Choose "dna", "rna" or "dnarna".')
    parser.add_argument('--subs_select', type=str, default="breast", help='Choose "all", "breast" or "lung".')
    parser.add_argument('--sams_select', type=str, default="bas",
                        help='Choose "all", "bas", "ont", "eot", "baseot", "basonteot".')
    parser.add_argument('--gene_only', type=str, default="yes",
                        help="If 'yes', simplify alteration names to gene names.")
    parser.add_argument('--col_side', type=str, default="response",
                        help="Choose between 'benefit', 'response', and 'timepoint'.")
    parser.add_argument('--run_tests', type=str, default="yes",
                        help="Choose 'no' to deactivate Fischer-Booschlo association tests.")
    parser.add_argument('--top_mode', type=str, default="alt",
                        help="Choose 'alt' to compute statistics only for alterations selected, 'all' for all.")
    parser.add_argument('--tests_alt', type=str, default="two-sided",
                        help="If 'greater', the test will assess whether the frequency of alterations is higher" + \
                             " in responder (resp. EOT) than in non-responder (resp. baseline). If 'less', the" + \
                             " other way around.")
    parser.add_argument('--output_pdf', type=str,  help='Path to output oncoplot.',
                default=f"{R_FOLDER}/breast/drivers_only/baseline_only_count_min_5_response.pdf")
    parser.add_argument('--output_tab', type=str,  help='Path to output table.',
                default=f"{R_FOLDER}/breast/drivers_only/baseline_only_count_min_5_response.xlsx")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)

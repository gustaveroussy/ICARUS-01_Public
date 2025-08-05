# -*- coding: utf-8 -*-
"""
@created: Jul 07 2022
@modified: Aug 29 2024
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Useful functions for drawing oncoplot-like figures.
"""

from functools import reduce
import numpy as np
import pandas as pd

# matplotlib
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# comut
import palettable
from comut import comut

# matplotlib
import matplotlib.cm as cm

# stats
from scipy.stats import mannwhitneyu, boschloo_exact
from statsmodels.sandbox.stats.multicomp import multipletests

def drop_duplicates_list(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def get_ordered_alterations_and_samples(df_cln, df_alt, sams, col_gen, col_alt, col_sub_id, col_sam_id, threshold_alt,
                                        threshold_alt_sams=None, complete_gene=False, contiguous_sub=False):
    """
    From a table of alteration data per sample, return a list of ordered samples for the
    columns of the oncoplot and a list or ordred alterations for the rows of the oncoplot.

    Parameters
    ----------
    df_alt: dataframe
        Must contain the columns `col_alt`, `col_gen`, `col_sub_id`, and `col_sam_id`.
    sams: list
        List of all samples. This is needed as some samples may be absent from the table `df_alt`.
    col_gen: str
        Name designating genes.
    col_alt: str
        Name designating alterations.
    col_sub_id: str
        Name designating the subject id. Useful if you want to position next to each other samples from the same
        subject. Otherwise may be identical to `col_sam_id`.
    col_sam_id: str
        Name designating the sample id.
    threshold_alt: double or int
        Value in [0,1) or int. Only alterations seen more frequently than this threshold will be retained.
    threshold_alt_sams: list
        List of samples to restricting to for the selection of recurrrent alterations
    complete_gene: bool, optional
        Set to True to list contiguously all alterations from the same in case the most frequent alteration of this gene
        passes the threshold `threshold_alt`.
    contiguous_sub: bool, optional
        Set to False to prevent reordering of samples so that samples coming from the same subject are contiguous.

    Returns
    -------
    alts_ordered, sams_ordered: tuple
        Tuple of 2 ordered lists.
    """
    # order the alterations from least frequent to most frequent
    # if complete_gene is set to True, the second, third, etc. most frequent alteration of the gene are
    # positioned after the most frequent alteration.
    if not complete_gene:
        if threshold_alt_sams is not None:
            df_alt_sams = df_alt.loc[df_alt[col_sam_id].isin(threshold_alt_sams)].copy()
        else:
            df_alt_sams = df_alt.copy()
        df_alt_cnt = df_alt_sams[[col_sam_id, col_alt]].drop_duplicates()[col_alt].value_counts()
        df_alt_cnt = df_alt_cnt.to_frame("Count").reset_index()
        df_alt_cnt = df_alt_cnt.rename(columns={"index": col_alt})
        df_alt_cnt = df_alt_cnt.sort_values(by=["Count", col_alt])
        df_alt_cnt = df_alt_cnt.set_index(col_alt)["Count"]
        if int(threshold_alt)==threshold_alt:
            alts_ordered = df_alt_cnt.loc[df_alt_cnt >= threshold_alt].index.tolist()
        else:
            alts_ordered = df_alt_cnt.loc[df_alt_cnt >= len(sams)*threshold_alt].index.tolist()
    else:
        # select only alterations seen more frequently than alt_threshold
        if threshold_alt_sams is not None:
            df_alt_sams = df_alt.loc[df_alt[col_sam_id].isin(threshold_alt_sams)].copy()
        else:
            df_alt_sams = df_alt.copy()
        df_alt_cnt = df_alt_sams[[col_sam_id, col_alt]].drop_duplicates()[col_alt].value_counts()
        df_alt_cnt = df_alt_cnt.to_frame("Count").reset_index()
        df_alt_cnt = df_alt_cnt.rename(columns={"index": col_alt})
        df_alt_cnt = df_alt_cnt.sort_values(by=["Count", col_alt])
        df_alt_cnt = df_alt_cnt.set_index(col_alt)["Count"]
        if int(threshold_alt)==threshold_alt:
            alts_ordered = df_alt_cnt.loc[df_alt_cnt >= threshold_alt].index.tolist()
        else:
            alts_ordered = df_alt_cnt.loc[df_alt_cnt >= len(sams)*threshold_alt].index.tolist()

        df_alt_sub = df_alt.loc[df_alt[col_alt].isin(alts_ordered)]
        alts_orderings = df_alt_sub[col_gen].value_counts().index[::-1].tolist()

        for alt_gen in alts_orderings:
            df_alt_gen = df_alt.loc[df_alt[col_gen]==alt_gen]
            alts_ordered_gen = df_alt_gen[col_alt].value_counts().index.tolist()[::-1]
            if len(alts_ordered_gen)>1:
                index_highest = alts_ordered.index(alts_ordered_gen[-1]) - len(alts_ordered)
                alts_ordered = [x for x in alts_ordered if x not in alts_ordered_gen]
                for i, alt in enumerate(alts_ordered_gen):
                    index_insert = index_highest + i + 1
                    if index_insert==0:
                        alts_ordered += [alt]
                    else:
                        alts_ordered.insert(index_insert, alt)
                    index_highest -= 1

    # samples are ordered using stratified lexicographical ordering
    cols_sort = [col_alt, col_sam_id, col_sub_id]
    cols_sort = drop_duplicates_list(cols_sort)
    df_alt_sam = df_alt[cols_sort].drop_duplicates()
    df_alt_sam["Mutated"] = "1"
    df_alt_binary = df_alt_sam.pivot(index=col_alt, columns=col_sam_id, values="Mutated")
    df_alt_binary = df_alt_binary.loc[alts_ordered[::-1],:].fillna("0")
    df_alt_string = df_alt_binary.apply(lambda x: "".join(x.tolist()), axis=0)
    sams_alts_ordered = df_alt_string.sort_values(ascending=False).index.tolist()
    sams_alts_ordered = sams_alts_ordered + [x for x in sams if x not in sams_alts_ordered]

    # samples with 0 detected alteration are appended at the end
    sams_no_alt = sorted(list(set(sams).difference(set(sams_alts_ordered))))
    sams_ordered = sams_alts_ordered + sams_no_alt

    # if 2 samples come from the same subject, there are place next to each other
    if contiguous_sub:
        sams_orderings = df_cln[col_sub_id].unique().tolist()

        for subject in sams_orderings:
            df_subject = df_cln.loc[df_cln[col_sub_id]==subject]
            sams_ordered_sub = sorted(df_subject[col_sam_id].unique().tolist())
            if len(sams_ordered_sub)>1:
                indexes_all =[sams_ordered.index(sam) for sam in sams_ordered_sub]
                index_highest = min(indexes_all) - len(sams_ordered)
                sams_ordered = [x for x in sams_ordered if x not in sams_ordered_sub]
                index_highest += len(sams_ordered_sub)
                for i, sam in enumerate(sams_ordered_sub):
                    index_insert = index_highest + i
                    if index_insert==0:
                        sams_ordered += [sam]
                    else:
                        sams_ordered.insert(index_insert, sam)
                    index_highest -= 1

    print("-%s alerations will be displayed" % len(alts_ordered))
    print("-%s samples will be displayed" % len(sams_ordered))

    return alts_ordered, sams_ordered


def add_responder_status(df_cln):
    """
    Add a `Responder` column that will serve to order samples and that will be displayed in the plots.

    Parameters
    ----------
    df_cln: dataframe
        Must contain the columns `Best_Overall_Response`.

    Returns
    -------
    df_cln: identical to input `df_cln` with an extra column.
    """
    mask_cr = (df_cln["Best_Overall_Response"]=="CR")
    mask_pr_c = (df_cln["Best_Overall_Response"]=="PR conf.")
    mask_pr_u = (df_cln["Best_Overall_Response"]=="PR unconf.")
    mask_sd = (df_cln["Best_Overall_Response"]=="SD")
    mask_pd = (df_cln["Best_Overall_Response"]=="PD")

    df_cln["Responder"] = "N/A"
    df_cln.loc[mask_cr, "Responder"] = "Yes"
    df_cln.loc[mask_pr_c, "Responder"] = "Yes"
    df_cln.loc[mask_pr_u, "Responder"] = "No"
    df_cln.loc[mask_sd, "Responder"] = "No"
    df_cln.loc[mask_pd, "Responder"] = "No"
    return df_cln


def get_tables(df_cln, df_alt, col_alt, col_alt_cat, col_sub_id, col_sam_id, alts_ordered, df_alt_top=None,
               df_alt_side=None, df_cln_side=None, categorical_dnv={}, col_side_barplot="Cohort",
               values_side_barplot_ignore=[]):
    """
    Prepare a dict of tables for the oncoplot.

    Parameters
    ----------
    df_cln: dataframe
        Must contain the columns `col_sub_id` and `col_sam_id`.
    df_alt: dataframe
        Must contain the columns `col_alt`, `col_gen`, `col_sub_id`, and `col_sam_id`.
    col_alt: str
        Name designating alterations.
    col_alt_cat: str
        Name designating alteration class. Used for the coloring of cells in the oncoplot.
    col_sub_id: str
        Name designating the subject id. Useful if you want to position next to each other samples from the same
        subject. Otherwise may be identical to `col_sam_id`.
    col_sam_id: str
        Name designating the sample id.
    alts_ordered: list
        List of ordered alterations.
    df_alt_top: dataframe
        Data for the top barplot. If not specified, default to all alterations in df_alt except for "VAF alterations".
    df_alt_side: dataframe
        Data for the side barplot. If not specified, default to alterations in the plot except for "VAF alterations".
    df_cln_side: dataframe
        Data for total of the side barplot. If not specified, default to `df_cln`.
    categorical_dnv: dict
        Keys are fixed tags defined in the code like 'her2', 'ben', etc. while values are the labels to appear on the
        plot.

    Returns
    -------
    dfs_data: dict
        Dict of tables for the oncoplot.
    """

    dfs_data = {}
    cols_comut = ["sample", "category", "value"]
    df_alt_no_vaf = df_alt.loc[~df_alt[col_alt_cat].apply(lambda x: "VAF" in x if type(x)==str else False)]
    df_alt_no_vaf_inplot = df_alt_no_vaf.loc[df_alt_no_vaf[col_alt].isin(alts_ordered)]

    # dataframe for indicator same patient
    df_cnt = df_cln.groupby(col_sub_id)[col_sam_id].size()
    subs_mult = df_cnt.loc[df_cnt > 1].index.tolist()
    df_ind = df_cln.loc[df_cln[col_sub_id].isin(subs_mult)][[col_sam_id, col_sub_id]]
    df_ind = df_ind.rename(columns={col_sam_id: "sample", col_sub_id: "group"})
    map_grp = {sub: i for i, sub in enumerate(subs_mult)}
    df_ind["group"] = df_ind["group"].map(map_grp)
    dfs_data["ind"] = df_ind

    # dataframe for center comut plot
    df_cen = df_alt.rename(columns={col_sam_id: "sample", col_alt: "category", col_alt_cat: "value"})
    df_cen = df_cen[cols_comut].drop_duplicates()
    dfs_data["cen"] = df_cen

    # dataframe for categorical data cohort
    df_coh = df_cln.rename(columns={col_sam_id: "sample", "Cohort": "value"})
    if "coh" in categorical_dnv:
        df_coh["category"] = categorical_dnv["coh"]["name"]
    else:
        df_coh["category"] = "Cohort"
    df_coh = df_coh[cols_comut].drop_duplicates()
    dfs_data["coh"] = df_coh

    for d, nv in categorical_dnv.items():
        df_cat = df_cln.rename(columns={col_sam_id: "sample", nv["var"]: "value"})
        df_cat["category"] = nv["name"]
        df_cat = df_cat[cols_comut].drop_duplicates()
        dfs_data[d] = df_cat

    # dataframe for top barplot
    if df_alt_top is None:
        df_alt_top = df_alt_no_vaf
    df_bur = df_alt_top.groupby(col_sam_id)[col_alt_cat].value_counts().unstack(level=-1)
    df_bur = df_bur.fillna(0).astype(int).reset_index()
    df_bur = df_bur.rename(columns={col_sam_id: "sample"})
    dfs_data["bur"] = df_bur

    # dataframe for side barplot
    if df_alt_side is None:
        df_alt_side = df_alt_no_vaf_inplot[[col_sam_id, col_alt, col_side_barplot]].drop_duplicates()
    if df_cln_side is None:
        df_cln_side = df_cln
    df_cnt = df_alt_side.groupby([col_alt, col_side_barplot]).size().unstack(level=-1)
    s_tot = df_cln_side[[col_sam_id, col_side_barplot]].drop_duplicates()[col_side_barplot].value_counts()
    df_frq = 100*df_cnt/s_tot
    df_frq = df_frq.fillna(0)
    df_frq.columns = df_frq.columns.tolist()
    df_cnt.columns = df_cnt.columns.tolist()
    df_frq = df_frq.reset_index().rename(columns={col_alt: "category"})
    df_cnt = df_cnt.reset_index().rename(columns={col_alt: "category"})

    # if some values should be ignored
    df_frq = df_frq[[x for x in df_frq if x not in values_side_barplot_ignore]]
    df_cnt = df_cnt[[x for x in df_cnt if x not in values_side_barplot_ignore]]

    dfs_data["frq"] = df_frq
    dfs_data["cnt"] = df_cnt

    # dataframe for stars
    df_sym = df_alt_no_vaf_inplot.rename(columns={col_sam_id: "sample", col_alt: "category", "Annotated": "value"})
    df_sym = df_sym[["sample", "category", "value"]]
    df_sym = df_sym.loc[df_sym["value"].apply(lambda x: "Yes" in str(x))].copy()
    df_sym["value"] = "Annotated OncoKB/CIViC"
    dfs_data["sym"] = df_sym

    return dfs_data


def get_mappings(t_vaf_labs, name_t_vaf_inc=None, borders=[], darkgrey_frq=False, col_side_barplot="Cohort",
                 labs_side_barplot=None, data_side_barplot=None):
    cmap_tab20 = cm.get_cmap("tab20")
    cmap_tab20b = cm.get_cmap("tab20b")
    cmap_tab10 = cm.get_cmap("tab10")
    cmap_green = cm.get_cmap("Greens")

    mappings = {}
    mappings["cen"] = {"Mutation": "#FFBC42",
                       "Indel": "#8F2D56",
                       "Amplification": "#D81159",
                       "HL Amp": "#D81159",
                       "ML Amp": "#ff4d6d",
                       "LL Amp": "#ff8fa3",
                       "LOH": "#48cae4",
                       "Deletion": "#0496FF",
                       "Hom. Deletion": "#0496FF",
                       "Hom Del": "#0496FF",
                       "Multiple": "#BC8034",
                       "Additional Mutation": {'facecolor':'none', 'edgecolor':'#BC8034', 'linewidth': 3},
                       "Additional Indel": {'facecolor':'none', 'edgecolor':'#362023', 'linewidth': 3},
                       "Additional Deletion": {'facecolor':'none', 'edgecolor':'#247BA0', 'linewidth': 3},
                       "Additional Amplification": {'facecolor':'none', 'edgecolor': '#B20D30', 'linewidth': 3}}

    if name_t_vaf_inc is not None:
        mappings["cen"][name_t_vaf_inc] = {'facecolor':'none', 'edgecolor':'#04724D', 'linewidth': 3}

    for i, t_vaf_lab in enumerate(t_vaf_labs):
        color_rgba = list(cmap_green(0.1 + 0.899 * i/(len(t_vaf_labs)-1)))
        color_rgba = [round(x*255) for x in color_rgba]
        color_hex = '#{:02x}{:02x}{:02x}'.format(*tuple(color_rgba))
        mappings["cen"][t_vaf_lab] = color_hex

    mappings["coh"] = {"BREAST-BAS": "#FFE381", "BREAST-EOT": "#CEC288", "TCGA-BRCA": "#9F9FED",
                       "LUNG-BAS": "#FFE381", "LUNG-EOT": "#CEC288"}
    mappings["tim"] = {"BAS": "#FFE381", "EOT": "#CEC288"}
    mappings["can"] = {"BREAST": "#FFC8DD", "LUNG": "#DBDB8D"}

    mappings["ben"] = {"Yes": "#0BC9CD", "No": "#FF9B85", "N/A": "lightgrey", "NE": "lightgrey"}
    mappings["res"] = {"CR": "#AAF683", "PR": "#60D394", "PR conf.": "#60D394",  "PR unconf.": "#dddf00",
                       "SD w/ benefit": "#FFD97D", "SD w/o benefit": "#FF9B85", "SD": "#FFDAB9", "PD": "#EE6055",
                       "N/A": "lightgrey", "NE": "lightgrey"}
    mappings["res_bin"] = {"Yes": "#60D394", "No": "#FF9B85"}
    mappings["mtc"] = {"Matched normal": "#D7BCE8", "Unmatched": "#8884FF"}
    mappings["bur"] = {k:v for k,v in mappings["cen"].items() if k not in borders}
    mappings["his"] = {"LUAD": "#bde0fe", "LUSC": "#8D5B4C", "LUAS": "#2F131E", "Lung": "lightgrey", "N/A": "lightgrey",
                       "ILC": "#7284A8", "IDC": "#7CDEDC", "NSQ": "#bde0fe", "SCC": "#8D5B4C",
                       "BRCA": "pink", "SPC": "#bde0fe", "MDLC": "#8D5B4C"}

    if labs_side_barplot is not None:
        mappings[data_side_barplot] = {}
        for i, lab_side_barplot in enumerate(labs_side_barplot):
            if len(labs_side_barplot) <= 10:
                color_rgba = list(cmap_tab10(0.1 + 0.899 * i/(len(labs_side_barplot)-1)))
            else:
                color_rgba = list(cmap_tab20(0.1 + 0.899 * i/(len(labs_side_barplot)-1)))
            color_rgba = [round(x*255) for x in color_rgba]
            color_hex = '#{:02x}{:02x}{:02x}'.format(*tuple(color_rgba))
            mappings[data_side_barplot][lab_side_barplot] = color_hex


    if darkgrey_frq:
        mappings["frq"] = {"Number of samples": "darkgrey"}
    else:
        if col_side_barplot=="Cohort":
            mappings["frq"] = mappings["coh"]
        elif col_side_barplot=="Clinical_Benefit":
            mappings["frq"] = mappings["ben"]
        elif col_side_barplot=="Responder":
            mappings["frq"] = mappings["res_bin"]
        elif col_side_barplot in ["Response_After_3_Months", "Best_Overall_Response"]:
            mappings["frq"] = mappings["res"]
        elif col_side_barplot=="Timepoint":
            mappings["frq"] = mappings["tim"]
        elif col_side_barplot=="Batch":
            mappings["frq"] = mappings[data_side_barplot]
        else:
            raise NotImplementedError("Please add colors for side barplot of %s" % col_side_barplot)

    mappings["sym"] = {"Annotated OncoKB/CIViC": ("black", "black")}

    return mappings


def add_vaf_for_mutations(df_alt, df_cln, col_sub_id, col_sam_id, col_alt, col_alt_det, col_alt_cat, col_t_vaf,
                          col_alt_cat_sec, mode="tex"):
    """
    If the VAF is represented in the oncoplot colors, then only one mutation may be shown for each combination of row
    and column. In this case, let us select the mutation with the highest VAF. Exception to this rule:
      - if there are 2 alterations or more and there is another paired sample,  prefer the mutation that is also seen
      in the paired sample

    Parameters
    ----------
    df_cln: dataframe
        Must contain the columns `col_sub_id` and `col_sam_id`.
    df_alt: dataframe
        Must contain the columns `col_alt`, `col_gen`, `col_sub_id`, and `col_sam_id`.
    col_sub_id: str
        Name designating the subject id. Useful if you want to position next to each other samples from the same
        subject. Otherwise may be identical to `col_sam_id`.
    col_sam_id: str
        Name designating the sample id.
    col_alt: str
        Name designating alterations.
    col_alt_det: str
        Name designating detailed alterations.
    col_alt_cat: str
        Name designating alteration class. Used for the coloring of cells in the oncoplot.
    col_t_vaf: str
        Name designating the VAF. Used for the coloring of cells in the oncoplot.
    col_alt_cat_sec: str
        Name designating secondary alteration class. Used for the coloring of cells in the oncoplot.
    mode: str
        Choose 'tex' or 'unicode'.

    Returns
    -------
    df_alt, t_vaf_labs: tuple dataframe, list
        Updated table with 1 additional row for each mutation/indel with the VAF category as a value for `col_alt_cat` and
        list of VAF labels.
    """
    df_cln_gby = df_cln.groupby([col_sub_id])["Cohort"].nunique()
    subjects_paired = df_cln_gby[df_cln_gby > 1].index.tolist()
    samples_paired = df_cln.loc[df_cln[col_sub_id].isin(subjects_paired)][col_sam_id].unique().tolist()
    mask_double = df_alt[[col_sam_id, col_alt, col_alt_cat]].duplicated(keep=False)
    mask_paired = df_alt[col_sam_id].isin(samples_paired)

    # mask a is for alterations repeated in the same sample for which another matched sample is also present
    mask_a = mask_double & mask_paired

    # mask b is for alterations repeated in the same sample for which no matched sample is also present
    mask_b = mask_double & ~mask_paired

    # mask c is for alterations not repeated in the same sample
    mask_c = ~mask_double
    df_alt_a = df_alt.loc[mask_a]
    df_alt_b = df_alt.loc[mask_b]
    df_alt_c = df_alt.loc[mask_c]

    # in case an alterations is repeated multiple times across multiple samples, select the one alteration with
    # the highest recurrence and, if ties left, with the highest VAF
    df_alt_a_cnt = df_alt_a.groupby([col_sub_id, col_alt, col_alt_det]).size().to_frame("Count_Alt").reset_index()
    df_alt_a_uni = df_alt_a.merge(df_alt_a_cnt, how="left", on=[col_sub_id, col_alt, col_alt_det])
    df_alt_a_uni = df_alt_a_uni.sort_values(by=["Count_Alt", col_alt, col_alt_det, "t_vaf"])
    df_alt_a_uni = df_alt_a_uni.drop_duplicates(subset=[col_sam_id, col_alt], keep="first")

    # in case an alterations is repeated multiple times across one sample, select the one alteration with
    # with the highest VAF
    df_alt_b_uni = df_alt_b.sort_values(by=[col_alt, "t_vaf"], ascending=False)
    df_alt_b_uni = df_alt_b_uni.drop_duplicates(subset=[col_sam_id, col_alt], keep="first")

    df_alt = pd.concat((df_alt_a_uni, df_alt_b_uni, df_alt_c))

    # create column col_alt_cat that will be the basis for the colors of the oncoplot
    # in order to allow a double annotation of VAF category and alteration identity for mutations,
    # add artificial rows where col_alt_cat will the rows
    t_vaf_bins = [0, 0.05-1e-5, 0.1, 0.2, 0.3, 0.4, 0.5, 1+ 1e-5]
    t_vaf_labs = []
    for i, t_vaf_bin in enumerate(t_vaf_bins[:-1]):
        if i == 0:
            if mode=="tex":
                t_vaf_labs.append("VAF $<$ %.2g" % (t_vaf_bins[i+1]))
            else:
                t_vaf_labs.append(u"VAF < %.2g" % (t_vaf_bins[i+1]))
        elif i < len(t_vaf_bins)-2:
            if mode=="tex":
                t_vaf_labs.append("%.2g $\\leq$ VAF $<$ %.2g" % (t_vaf_bins[i], t_vaf_bins[i+1]))
            else:
                t_vaf_labs.append(u"%.2g ≤ VAF < %.2g" % (t_vaf_bins[i], t_vaf_bins[i+1]))
        else:
            if mode=="tex":
                t_vaf_labs.append("%.2g $\\leq$ VAF $\\leq$ 1" % (t_vaf_bins[i]))
            else:
                t_vaf_labs.append("%.2g ≤ VAF ≤ 1" % (t_vaf_bins[i]))

    mask_mut = df_alt[col_alt_cat].isin(["Mutation", "Indel"])
    mask_cat_sec = ~df_alt[col_alt_cat_sec].isnull()

    # for alterations that are mutations only (not combined with alterations of another nature)
    # add a secondary row just for the VAF
    df_alt_mut = df_alt.loc[mask_mut & ~mask_cat_sec].copy()
    df_alt_mut[col_alt_cat] = pd.cut(df_alt_mut[col_t_vaf], bins=t_vaf_bins, labels=t_vaf_labs)

    # for alterations that combine events of different natures, add a secondary row for the secondary nature
    df_alt_cat_sec = df_alt.loc[mask_cat_sec].copy()
    df_alt_cat_sec[col_alt_cat] = df_alt_cat_sec[col_alt_cat_sec]

    df_alt = pd.concat((df_alt_cat_sec, df_alt_mut, df_alt), axis=0)

    return df_alt, t_vaf_labs


def convert_num_to_str(x):
    try:
        if int(x)==x:
            y = "%d" % x
        else:
            raise ValueError
    except:
        try:
            y = "%.1f" % x
            if y=="nan":
                y = x
        except:
            y = x

    return y


def draw_comut_plot(dfs_data, mappings, sams_ordered, alts_ordered, t_vaf_labs, categorical_data=[],
                    categorical_names=[], borders=None, width_left=None, width_middle=None, shadow_width_left=None,
                    height_middle=None, height_top=None, hspace=0.01, wspace=0.05, xtick_fontsize=8, ytick_fontsize=12,
                    label_bar_top="Drivers", label_bar_top_fontsize=12, stacked_bar_top=True, stacked_bar_side=True,
                    label_bar_side="Number of samples", label_bar_side_fontsize=12, digits_bar_size_fontsize=None,
                    mode_bar_side="Number", sample_indicators=False, ncol_legend=1, bbox_to_anchor_legend=(1,1),
                    ignored_plots=["Alteration Type", "Same patient bis"], ignored_values=[], gap_between_groups=0.2,
                    gap_within_groups=0.05, labels_orders={}, show_bar_side=True, show_symbols=False,
                    show_legend_barplot=True, title_legend_barplot=None, unit_width=0.15, unit_height=0.20):
    """
    Instantiate a `CoMut` object and populate with data. Graphical parameters are adjusted based on user parameters.


    Parameters
    ----------
    dfs_data: dict
        Dict of dataframes with data for the figure.
    mappings: dict
        Dict of colors. Keys should match keys of `dfs_data`
    sams_ordered: list
        Ordered samples.
    alts_ordered: list
        Ordered alterations.
    t_vaf_labs: list
        Names of the VAF classes.
    categorical_data: list
        List of keys of `dfs_data` holding categorical data to be displayed at the top.
    categorical_names: list
        List of the labels that will be displayed for the categorical data.
    ncol_legend: int
        Number of columns for legend.


    Returns
    ----------
    comut_alt: a `CoMut` object
        Populated with data and parameters set for the rendering the figure.
    """

    comut_alt = comut.CoMut()
    comut_alt.samples = sams_ordered

    priority = t_vaf_labs[::-1] + ["Mutation", "Amplification", "HL Amp", "Deletion", "Hom. Deletion", "Hom Del",
                                   "Indel", "ML Amp", "LL Amp", "LOH"]
    value_order = t_vaf_labs[::-1] + ["Mutation", "Amplification", "HL Amp", "Deletion", "Hom. Deletion",
                                      "Hom Del", "Indel", "ML Amp", "LL Amp", "LOH"]

    if sample_indicators:
        # add indicators first, since they will be at the bottom
        indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 5}

        comut_alt.add_sample_indicators(data=dfs_data["ind"], name='Same patient', plot_kwargs=indicator_kwargs,
                                        xtick_show=True, xtick_fontdict={"fontsize": xtick_fontsize})

        comut_alt.add_categorical_data(data=dfs_data["cen"], name='Alterations', category_order=alts_ordered,
                                       borders=borders, priority=priority, value_order=value_order,
                                       mapping=mappings["cen"], xtick_show=False,
                                       xtick_fontdict={"fontsize": xtick_fontsize}, ytick_style='italic',
                                       ytick_fontdict={"fontsize": ytick_fontsize})

        comut_alt.add_sample_indicators(data=dfs_data["ind"], name='Same patient bis', plot_kwargs=indicator_kwargs,
                                        xtick_show=True, xtick_fontdict={"fontsize": 8})
    else:
        comut_alt.add_categorical_data(data=dfs_data["cen"], name='Alterations', category_order=alts_ordered,
                                       borders=borders, priority=priority, value_order=value_order,
                                       mapping=mappings["cen"], xtick_show=True,
                                       xtick_fontdict={"fontsize": xtick_fontsize}, ytick_style='italic',
                                       ytick_fontdict={"fontsize": ytick_fontsize})

    for cat_data, cat_name in zip(categorical_data, categorical_names):
        comut_alt.add_categorical_data(data=dfs_data[cat_data], name=cat_name,
                                       mapping=mappings[cat_data], xtick_show=False,
                                       ytick_style='normal', ytick_fontdict={"fontsize": ytick_fontsize})

    comut_alt.add_bar_data(data=dfs_data["bur"], name='Alteration Type',
                           mapping=mappings["bur"], ytick_fontdict={"fontsize": 12}, stacked=stacked_bar_top,
                           ylabel=label_bar_top, ylabel_fontsize=label_bar_top_fontsize, ylabel_rotation=90)

    if show_bar_side:
        if mode_bar_side=="Number":
            data = dfs_data["cnt"]
        elif mode_bar_side=="Percent":
            data = dfs_data["frq"]

        comut_alt.add_side_bar_data(data=data, paired_name='Alterations', name="Number of samples",
                                    position="left", mapping=mappings["frq"], xtick_fontdict={"fontsize": 12},
                                    stacked=stacked_bar_side, xlabel=label_bar_side,
                                    xlabel_fontsize=label_bar_side_fontsize, xlabel_rotation=0,
                                    gap_between_groups=gap_between_groups, gap_within_groups=gap_within_groups)

    if show_symbols:
        comut_alt.add_scatter_data(data=dfs_data["sym"], paired_name='Alterations', name='Annotations',
                                   mapping=mappings["sym"], marker="*", markersize=8)


    # widths
    if show_bar_side:
        if width_left is None:
            width_left = 4
        if width_middle is None:
            width_middle = unit_width*len(sams_ordered)
        if shadow_width_left is None:
            shadow_width_left = 2.5

        total_width = width_left + width_middle + shadow_width_left
        r_width_left = width_left/total_width
        r_width_middle = width_middle/total_width
        r_shadow_width_left = shadow_width_left/total_width
        widths = [r_width_left, r_width_middle]
    else:
        if width_middle is None:
            width_middle = unit_width*len(sams_ordered)
        total_width = width_middle
        r_shadow_width_left=None
        widths = None

    # heights
    if height_middle is None:
        height_middle = unit_height*len(alts_ordered)
    if height_top is None:
        height_top = max(1.75*unit_height*height_middle, 2.5)

    total_height = height_top + height_middle
    heights = {'Alteration Type': height_top}

    # render the plot
    comut_alt.plot_comut(x_padding=0.04,
                         y_padding=0.04,
                         tri_padding=0.03,
                         figsize=(total_width,total_height),
                         hspace=hspace,
                         wspace=wspace,
                         heights=heights,
                         widths=widths,
                         shadow_width_left=r_shadow_width_left)

    # configure legends
    handles_more = []
    labels_more = []
    titles_more = []

    # legend for VAF colors
    labels = t_vaf_labs
    colors = [mappings["cen"][lab] for lab in labels]
    handles = [mpatches.Patch(color=color) for color in colors]

    handles_more.append(handles)
    labels_more.append(labels)
    titles_more.append("Alteration VAF")

    # legend for cohorts
    df_frq = dfs_data["frq"]
    n_cohorts = df_frq.shape[1]-1
    if n_cohorts > 1 and show_legend_barplot:
        labels = df_frq.columns[1:].tolist()
        colors = [mappings["frq"][lab] for lab in labels]
        handles = [mpatches.Patch(color=color) for color in colors]

        handles_more.append(handles)
        labels_more.append(labels)
        titles_more.append(title_legend_barplot)

    if show_symbols:
        handles = []
        labels = []
        for label, color in mappings["sym"].items():
            if len(label.split("&"))==1:
                handles.append(mlines.Line2D([], [], color=color[0], marker='*', linestyle='None', markersize=10))
                labels.append(label)

        handles_more.append(handles)
        labels_more.append(labels)
        titles_more.append("Annotation")

    comut_alt.add_unified_legend(ncol=ncol_legend, axis_name="Alterations", ignored_plots=ignored_plots,
                                 ignored_values=ignored_values, handles_more=handles_more, labels_more=labels_more,
                                 titles_more=titles_more, bbox_to_anchor=bbox_to_anchor_legend,
                                 labels_orders=labels_orders)

    # last adjustments

    # show spines central comut
    axis_name = "Alterations"
    for loc in ['bottom']:
        comut_alt.axes[axis_name].spines[loc].set_visible(True)

    # labels barplot
    axis_name = "Alteration Type"
    text_a = comut_alt.axes[axis_name].get_yticklabels()[0].get_text()
    text_b = comut_alt.axes[axis_name].get_yticklabels()[-1].get_text()
    text_a = "%d" % int(round(float(text_a)))
    text_b = "%d" % int(round(float(text_b)))
    comut_alt.axes[axis_name].set_yticklabels([text_a, text_b])

    # calculate the percentage of samples with that gene mutated, rounding and adding a percent sign
    if show_bar_side:
        axis_name = "Number of samples"
        df_frq = dfs_data["frq"]
        df_frq = df_frq.set_index("category").loc[alts_ordered].reset_index()
        n_bars = df_frq.shape[1]-1
        height = (1-gap_between_groups-(n_bars-1)*gap_within_groups)/n_bars

        yticks = []
        yticklabels = []

        for i,col in enumerate(df_frq.columns[1:]):
            y_range = -0.5 + gap_between_groups/2 + (i+0.5) * height + i * gap_within_groups + \
                    np.arange(0.5, len(alts_ordered))
            yticks += y_range.tolist()
            percentages = df_frq[col].apply(convert_num_to_str)
            yticklabels += percentages.tolist()

        # set location of yticks
        comut_alt.axes[axis_name].set_yticks(yticks)

        # set labels of yticks
        comut_alt.axes[axis_name].set_yticklabels(yticklabels, ha="right")

        # set label size
        if digits_bar_size_fontsize is not None:
            labelsize = digits_bar_size_fontsize
        else:
            if n_bars==1:
                labelsize = 12
            elif n_bars==2:
                labelsize = 8
            elif n_bars >= 3:
                labelsize = 5

        # move the ytick labels inside the bar graph
        comut_alt.axes[axis_name].tick_params(axis='y', pad=0, length=0, labelsize=labelsize)

        # Make y axis visible (by default it is not visible)
        comut_alt.axes[axis_name].get_yaxis().set_visible(True)

        # move y axis ticks to the right
        comut_alt.axes[axis_name].yaxis.tick_right()

        # x axis ticks
        if mode_bar_side=="Number":
            x_max = dfs_data["cnt"].set_index("category").max().max()
        elif mode_bar_side=="Percent":
            x_max = dfs_data["frq"].set_index("category").max().max()
        x_max_str = convert_num_to_str(x_max)
        comut_alt.axes[axis_name].set_xticks([0, x_max])
        comut_alt.axes[axis_name].set_xticklabels([0, x_max_str], fontdict={"fontsize": 12}, style="normal", rotation=0)

    # move vertical adjustment of ytick labels of burden plot
    axis_name = "Alteration Type"

    # set labels of yticks
    yticklabels = comut_alt.axes[axis_name].get_yticklabels()
    comut_alt.axes[axis_name].set_yticklabels(yticklabels, va="bottom")

    return comut_alt


def run_fisher_boschloo_tests(df_cln, df_alt, col_side_barplot, tests_alt, sams, col_gen, col_alt, col_sub_id,
                              col_sam_id, threshold_alt_tests=0.04, alts_ordered=None):
    """
    Runs multiple Fischer-Boschloo tests to assess the association between the binary values depicted on the left side
    barplot of the comut plot and the presence or absence of alterations among the list of alterations shown in the
    plot.
    The p-values are corrected for multiple testing using the Benjamini-Hochberg procedure.
    """

    # get table of counts per category
    df_cnt = df_alt[[col_sam_id, col_alt, col_side_barplot]].drop_duplicates()
    df_cnt = df_cnt.groupby([col_alt, col_side_barplot]).size().unstack(level=-1).fillna(0)
    margins = df_cln[col_side_barplot].value_counts()

    # get ordered alterations
    if alts_ordered is None:
        df_alt_cnt = df_alt_tab[[col_sam_id, col_alt]].drop_duplicates()[col_alt].value_counts().to_frame("Count")
        df_alt_cnt = df_alt_cnt.reset_index(drop=False)
        df_alt_cnt = df_alt_cnt.sort_values(by=["Count", "index"], ascending=False)
        alts_ordered = df_alt_cnt.index.tolist()[::-1]

    df_alt_tests = pd.DataFrame(columns=["Alteration", "P_Value"])
    for alt in alts_ordered[::-1]:
        if col_side_barplot=="Clinical_Benefit":
            val_a = "Yes"
            val_b = "No"
        elif col_side_barplot=="Responder":
            val_a = "Yes"
            val_b = "No"
        elif col_side_barplot=="Timepoint":
            val_a = "EOT"
            val_b = "BAS"

        if val_a not in df_cnt:
            cnt_a = 0
        else:
            cnt_a = df_cnt.loc[alt, val_a]

        if val_b not in df_cnt:
            cnt_b = 0
        else:
            cnt_b = df_cnt.loc[alt, val_b]

        table = np.array([[cnt_a, cnt_b],
                            [margins[val_a]-cnt_a, margins[val_b]-cnt_b]])
        out_test = boschloo_exact(table=table, alternative=tests_alt, n=100)
        pvalue = out_test.pvalue

        df_alt_tests = pd.concat((df_alt_tests, pd.DataFrame({"Alteration": [alt], "P_Value": [pvalue]})))

    # adjust pvalues
    df_alt_tests["P_Value_Adj"] = multipletests(pvals=df_alt_tests["P_Value"], alpha=0.1, method="fdr_bh")[1]

    return df_alt_tests


def select_alt_from_gene_lists(df_alt, col_gen, gene_list_path, gene_list_names):
    """
    Select only alterations affecting genes from the gene lists specified.
    """

    df_gene_list = pd.read_excel(gene_list_path)
    gene_list_union = []
    for gene_list_name in gene_list_names:
        gene_list = df_gene_list[gene_list_name].dropna().drop_duplicates().tolist()
        gene_list = sorted(list(set([e.strip() for g in gene_list for e in g.split(",")])))
        print("-INFO: selecting %d genes from gene list %s" % (len(gene_list), gene_list_name))
        gene_list_union.append(gene_list)

    gene_list_union = sorted(reduce(lambda l1,l2: list(set(l1).union(set(l2))), gene_list_union))
    mask_gene_list = df_alt[col_gen].isin(gene_list_union)
    print("-INFO: selected %d/%d alterations from %d genes" % \
          (sum(mask_gene_list), len(mask_gene_list), len(gene_list_union)))

    return df_alt.loc[mask_gene_list,:].copy()


def select_inactivating_alt(df_alt, col_alt_cat):
    """
    Select only inactivating alterations.
    """

    mut_inactivating = ["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Splice_Region",
                        "Nonstop_Mutation", "Stop_Codon_Del"]
    mask_mut = (df_alt[col_alt_cat].isin(["Mutation", "Indel"])) & \
            (df_alt["Variant_Classification"].isin(mut_inactivating))

    cna_inactivating = ["Deletion", "LOH", "cnLOH"]
    mask_cna = df_alt[col_alt_cat].isin(cna_inactivating)

    # for annotated events, rescue predicted LOF
    mask_okb = df_alt["Oncokb_Mutation_Effect"].isin(["Loss-of-function", "Likely Loss-of-function"])

    # merge
    mask_inactivating = mask_mut | mask_cna | mask_okb
    print("-INFO: selected %d/%d inactivating alterations" % (sum(mask_inactivating), len(mask_inactivating)))

    return df_alt.loc[mask_inactivating,:].copy()

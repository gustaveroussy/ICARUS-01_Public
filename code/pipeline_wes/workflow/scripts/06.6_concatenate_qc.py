# -*- coding: utf-8 -*-
"""
@created: Jun 07 2023
@modified: Jun 07 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Prepare a table compiling QC metrics on the analysis of WES files.
"""

import argparse
import gzip
import io
import os
import numpy as np
import pandas as pd
import re
import sys
sys.path.append("functions")
import zipfile

# functions ============================================================================================================

def read_header(path, prefix):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    else:
        with open(path, "r") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    return header


def read_and_process_facets_vcf(filepath):
    basename = os.path.basename(filepath)
    tsample = basename.split("_vs_")[0]
    nsample = basename.split("_vs_")[1].split(".vcf")[0]
    header = read_header(filepath, prefix="##")

    cols_old = ["purity", "ploidy", "emflags"]
    cols_new = ["Purity", "Ploidy", "EM_Flags"]
    values = {}

    for line in header:
        for col in cols_old:
            pattern = "##%s=" % col
            if line.startswith(pattern):
                try:
                    if col!="emflags":
                        value = float(line.split(pattern)[1])
                    else:
                        value = line.split(pattern)[1].strip()
                except:
                    value = np.nan

                values[col] = [value]

    df = pd.DataFrame(values)
    df = df.rename(columns={o:n for o,n in zip(cols_old, cols_new)})

    df.insert(0, "Tumor_Sample_Id", tsample)
    df.insert(1, "Normal_Sample_Id", nsample)
    return df



def convert_num_to_str(x, digits=2):
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
            if y.startswith("nan"):
                y = x
        except:
            y = x
    return y


def convert_num_to_str_with_unit(x, unit=1e6, digits=2):
    unit_to_suffix = {0.01: "%%", 1: "%%", 1e3: "K", 1e6: "M", 1e9: "B"}

    try:
        suffix = unit_to_suffix[unit]
    except:
        raise ValueError("Choose one of the following units: %s" % "-".join(list(unit_to_suffix.keys())))

    try:
        y = float(x)/unit
        if digits is not None:
            y = ("%%.%df %s" % (digits, suffix)) % y
        else:
            y = "%f %s" % (y, suffix)
        if y.startswith("nan"):
            y = x
    except:
        try:
            y = str(x)
            if y.startswith("nan"):
                y = x
        except:
            y = x
    return y


def merge_sample_on_tnp(df_sam, df_tnp, col_sam, col_dna_t, col_dna_n):
    cols_to_merge = [x for x in df_sam if x!=col_sam]

    df_sam_dna_t = df_sam.rename(columns={**{col_sam: col_dna_t}, **{x: "%s_T" % x for x in cols_to_merge}})
    df_tnp[col_dna_t] = df_tnp[col_dna_t].fillna("NA")
    df_tnp = df_tnp.merge(df_sam_dna_t, how="left", on=col_dna_t)
    df_tnp[col_dna_t] = df_tnp[col_dna_t].replace({"NA": np.nan})

    if col_dna_n in df_tnp and df_tnp[col_dna_n].isnull().mean() < 1:
        df_sam_dna_n = df_sam.rename(columns={**{col_sam: col_dna_n}, **{x: "%s_N" % x for x in cols_to_merge}})
        df_tnp[col_dna_n] = df_tnp[col_dna_n].fillna("NA")
        df_tnp = df_tnp.merge(df_sam_dna_n, how="left", on=col_dna_n)
        df_tnp[col_dna_n] = df_tnp[col_dna_n].replace({"NA": np.nan})


    return df_tnp


def add_purity_ploidy(df_tnp, somatic_facets):
    if len(somatic_facets) > 0:
        print("-adding purity and ploidy from %d files..." % len(somatic_facets), end="")
        dfs_ppy = []
        for file_facets in somatic_facets:
            dfs_ppy.append(read_and_process_facets_vcf(file_facets))
        df_ppy = pd.concat(dfs_ppy)
        df_ppy["DNA_P"] = df_ppy[["Tumor_Sample_Id", "Normal_Sample_Id"]].fillna("NA").apply("_vs_".join, axis=1)
        cols_ppy = ["DNA_P", "Purity", "Ploidy", "EM_Flags"]
        df_tnp["WES_somatic_cna"] = 0
        df_tnp.loc[df_tnp["DNA_P"].isin(df_ppy["DNA_P"]), "WES_somatic_cna"] = 1
        df_tnp = df_tnp.merge(df_ppy[cols_ppy], how="left", on="DNA_P")
        print("done!")
    return df_tnp


def add_mutation_statistics(df_tnp, somatic_maf, somatic_maf_agg):
    if len(somatic_maf) > 0 or somatic_maf_agg is not None:
        if len(somatic_maf) > 0:
            print("-adding mutation statistics from %d files..." % len(somatic_maf), end="")
            tnp_ids = []
            dfs_maf = []
            for file_maf in somatic_maf:
                tnp_ids.append(file_maf.split(".maf")[0])
                skiprows = len(read_header(path=file_maf, prefix="##"))
                dfs_maf.append(pd.read_table(file_maf, sep="\t", skiprows=skiprows, low_memory=False))
            df_maf = pd.concat(dfs_maf)
            df_tnp["WES_somatic_maf"] = 0
            df_tnp.loc[df_tnp["DNA_P"].isin(tnp_ids), "WES_somatic_maf"] = 1

        elif somatic_maf_agg is not None:
            print("-adding mutation statistics from 1 aggregated file...", end="")
            skiprows = len(read_header(path=somatic_maf_agg, prefix="##"))
            df_maf = pd.read_table(somatic_maf_agg, sep="\t", skiprows=skiprows, low_memory=False)

        col_t, col_n = "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"
        df_maf["DNA_P"] = df_maf[[col_t, col_n]].fillna("NA").apply("_vs_".join, axis=1)

        # number of mutations
        df_maf_nmut = df_maf.groupby("DNA_P").size().to_frame("N_Mutations").reset_index()

        # vaf quantiles
        df_maf["t_vaf"] = df_maf["t_alt_count"]/df_maf["t_depth"]
        quantiles = [0, 0.25, 0.5, 0.75, 1]
        cols_vafq = ["VAF_quantile_%d%%" % (q*100) for q in quantiles]
        df_maf_vafq = df_maf.groupby("DNA_P")["t_vaf"].apply(lambda x: np.quantile(x, quantiles)).to_frame("VAF")
        df_maf_vafq = df_maf_vafq.reset_index()
        df_maf_vafq[cols_vafq] = df_maf_vafq["VAF"].apply(pd.Series)
        del df_maf_vafq["VAF"]

        # add to main table
        cols_cur = df_tnp.columns.tolist()
        df_tnp = df_tnp.merge(df_maf_nmut, how="left", on="DNA_P")
        df_tnp = df_tnp.merge(df_maf_vafq, how="left", on="DNA_P")

        # fill zeroes to samples processed but with 0 mutations
        if "WES_somatic_maf" in df_tnp:
            cols_new = list(set(df_tnp.columns).difference(set(cols_cur)))
            mask_wes = df_tnp["WES_somatic_maf"]==1
            df_tnp.loc[mask_wes, cols_new] = df_tnp.loc[mask_wes, cols_new].fillna(0)
        print("done!")

    return df_tnp


def add_qc_fastqc(df_tnp, qc_fastqc, orientation="R1"):
    if len(qc_fastqc)>0:
        print("-adding FASTQC statistics %s from %d files..." % (orientation, len(qc_fastqc)), end="")
        col_sam = "Sample_Id"
        col_fastqc_reads = "FASTQC_Reads"
        df_qc = pd.DataFrame(columns=[col_sam, col_fastqc_reads])
        for file_qc in qc_fastqc:
            file_qc_base = re.sub(".zip$", "", os.path.basename(file_qc))
            sample_id = file_qc_base.split("_%s" % orientation)[0]
            with zipfile.ZipFile(file_qc) as zf:
                with io.TextIOWrapper(zf.open(os.path.join(file_qc_base, "fastqc_data.txt"))) as f:
                    lines = [l.strip() for l in f.readlines()]
                    lines_header = lines[:20]
                    for line in lines_header:
                        if line.startswith("Total Sequences"):
                            fastqc_reads = int(line.split("\t")[1])
                    row_qc = {col_sam: [sample_id], col_fastqc_reads: [fastqc_reads]}
                    df_qc = pd.concat((df_qc, pd.DataFrame(row_qc)))

        # rename columns according to orientation
        df_qc = df_qc.rename(columns={col_fastqc_reads: "%s_%s" % (col_fastqc_reads, orientation)})

        # add to main table
        df_tnp = merge_sample_on_tnp(df_qc, df_tnp, col_sam=col_sam, col_dna_t="DNA_T", col_dna_n="DNA_N")
        print("done!")

    return df_tnp


def add_qc_ncm(df_tnp, qc_ncm):
    if len(qc_ncm) > 0:
        print("-adding NGSCheckMate metrics from %d files..." % len(qc_ncm), end="")
        col_ncm_status = "NCM_Status"
        col_ncm_correlation = "NCM_Correlation"
        df_qc = pd.DataFrame(columns=["DNA_P", col_ncm_status, col_ncm_correlation])
        for file_qc in qc_ncm:
            tnp  = os.path.basename(file_qc).split("_all.txt")[0]
            with open(file_qc) as f:
                lines = [x.strip() for x in f.readlines()]
            qc_ncm_status = lines[0].split("\t")[1]
            qc_ncm_correlation = float(lines[0].split("\t")[3])
            row_qc = {"DNA_P": [tnp], col_ncm_status: [qc_ncm_status], col_ncm_correlation: [qc_ncm_correlation]}
            df_qc = pd.concat((df_qc, pd.DataFrame(row_qc)))

        # add to main table
        df_tnp = df_tnp.merge(df_qc, how="left", on="DNA_P")
        print("done!")

    return df_tnp


def add_qc_flagstat(df_tnp, qc_flagstat):
    if len(qc_flagstat) > 0:
        print("-adding flagstat metrics from %d files..." % len(qc_flagstat), end="")
        col_sam = "Sample_Id"
        col_flagstat_reads = "Flagstat_Passed_Reads"
        df_qc = pd.DataFrame(columns=["Sample_Id", col_flagstat_reads])

        for file_qc in qc_flagstat:
            with open(file_qc, "r") as f:
                lines = [x.strip() for x in f.readlines()]
                line = lines[0]
            sample_id = os.path.basename(file_qc).split("_flagstat")[0]
            qc_passed_reads = int(lines[0].split("+")[0].strip())
            row_qc = {col_sam: [sample_id], col_flagstat_reads: [qc_passed_reads]}
            df_qc = pd.concat((df_qc, pd.DataFrame(row_qc)))

        # add to main table
        df_tnp = merge_sample_on_tnp(df_sam=df_qc, df_tnp=df_tnp, col_sam=col_sam, col_dna_t="DNA_T", col_dna_n="DNA_N")
        print("done!")

    return df_tnp


def add_qc_hsmetrics(df_tnp, qc_hsmetrics):
    if len(qc_hsmetrics) > 0:
        print("-adding CollectHsMetrics metrics from %d files..." % len(qc_hsmetrics), end="")
        col_sam = "Sample_Id"
        df_qcs = []
        for file_qc in qc_hsmetrics:
            sample_id = os.path.basename(file_qc).split("_hs_metrics.tsv")[0]
            with open(file_qc, "r") as f:
                lines = [x.strip() for x in f.readlines()]
            header = lines[6].split("\t")
            values = lines[7].split("\t")
            df_qc = pd.DataFrame({h: [v] for h,v in zip(header, values)})
            df_qc[col_sam] = sample_id
            df_qcs.append(df_qc)
        df_qc = pd.concat(df_qcs)

        # select qc columns
        cols_qc = [col_sam, "TOTAL_READS", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE", "ZERO_CVG_TARGETS_PCT",
                   "PCT_TARGET_BASES_1X", "PCT_TARGET_BASES_2X", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X",
                   "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_40X", "PCT_TARGET_BASES_50X",  "PCT_TARGET_BASES_100X",
                   "HET_SNP_SENSITIVITY"]
        df_qc = df_qc[cols_qc]

        for col_qc in cols_qc[1:]:
            df_qc[col_qc] = df_qc[col_qc].astype(float)

        # add to main table
        df_tnp = merge_sample_on_tnp(df_sam=df_qc, df_tnp=df_tnp, col_sam=col_sam, col_dna_t="DNA_T", col_dna_n="DNA_N")
        print("done!")

    return df_tnp


def add_auto_qc(df_tnp, args):
    print("-adding 'Comment_QC_Auto' column...", end="")

    # add auto QC
    masks = []
    comments = []

    if "NCM_Status" in df_tnp:
        masks.append(df_tnp["NCM_Status"]=="unmatched")
        comments.append("Mismatch tumor/normal")

    if "MEDIAN_TARGET_COVERAGE_T" in df_tnp:
        masks.append(df_tnp["MEDIAN_TARGET_COVERAGE_T"].astype(float) < args.min_cov_t)
        comments.append("Median target cov in tumor < %s" % args.min_cov_t)

    if "MEDIAN_TARGET_COVERAGE_N" in df_tnp:
        masks.append(df_tnp["MEDIAN_TARGET_COVERAGE_N"].astype(float) < args.min_cov_n)
        comments.append("Median target cov in normal < %s" % args.min_cov_n)

    if "PCT_TARGET_BASES_10X_T" in df_tnp:
        masks.append(df_tnp["PCT_TARGET_BASES_10X_T"].astype(float) < args.min_10x_pct_t/100)
        comments.append("Target cov 10X in tumor < %s%%" % args.min_10x_pct_t)

    if "PCT_TARGET_BASES_10X_N" in df_tnp:
        masks.append(df_tnp["PCT_TARGET_BASES_10X_N"].astype(float) < args.min_10x_pct_n/100)
        comments.append("Target cov 10X in normal < %s%%" % args.min_10x_pct_n)

    if len(masks) > 0:
        df_tnp["Comment_QC_Auto"] = np.nan
        for mask, comment in zip(masks, comments):
            mask_null = df_tnp["Comment_QC_Auto"].isnull()
            df_tnp.loc[mask & mask_null, "Comment_QC_Auto"] = comment
            df_tnp.loc[mask & ~mask_null, "Comment_QC_Auto"] += " & %s" % comment
    print("done!")

    return df_tnp


def main(args):
    # load tumor normal pairs
    df_tnp = pd.read_table(args.tnp)

    # qc metrics from variant calling ==================================================================================

    # add estimated purity and ploidy, if any
    df_tnp = add_purity_ploidy(df_tnp, args.somatic_facets)

    # add number mutations and vaf quantiles, if any
    df_tnp = add_mutation_statistics(df_tnp, args.somatic_maf, args.somatic_maf_agg)

    # qc metrics pre-alignment =========================================================================================

    df_tnp = add_qc_fastqc(df_tnp, args.qc_fastqc_r1, orientation="R1")
    df_tnp = add_qc_fastqc(df_tnp, args.qc_fastqc_r2, orientation="R2")

    # qc metrics post-alignment ========================================================================================

    # add NGSCheckMate qc, if any
    df_tnp = add_qc_ncm(df_tnp, args.qc_ncm)

    # add number of qc-passed reads from samtools flagstat, if any
    df_tnp = add_qc_flagstat(df_tnp, args.qc_flagstat)

    # add many metrics computed by CollectHsMetrics (GATK tool)
    df_tnp = add_qc_hsmetrics(df_tnp, args.qc_hsmetrics)

    # last steps =======================================================================================================

    # add column "Comment_QC_Auto" with predefined thresholds
    df_tnp = add_auto_qc(df_tnp, args)

    if not args.not_human:
        # format reads column with million (M) unit
        cols_reads = [x for x in df_tnp if "reads" in x.lower()]
        for col_read in cols_reads:
            df_tnp[col_read] = df_tnp[col_read].apply(convert_num_to_str_with_unit, unit=1e6, digits=1)

        # format percent column with percent (%) unit
        cols_pct = [x for x in df_tnp if "PCT" in x or "SNP_SENSITIVITY" in x  or "quantile" in x]
        for col_pct in cols_pct:
            df_tnp[col_pct] = df_tnp[col_pct].apply(convert_num_to_str_with_unit, unit=0.01, digits=2)

    # format other int column
    cols_int = [x for x in df_tnp if x.startswith("N_") or "median" in x.lower()]
    for col_int in cols_int:
        df_tnp[col_int] = df_tnp[col_int].apply(convert_num_to_str)

    # save table of summary analysis
    if args.output.endswith(".xlsx"):
        df_tnp.to_excel(args.output, index=False, float_format="%.2f")
    elif args.output.endswith(".tsv"):
        df_tnp.to_csv(args.output, index=False, sep="\t", float_format="%.2f")
    else:
        df_tnp.to_csv(args.output, index=False, float_format="%.2f")
    print("-file saved at %s" % args.output)


if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Prepare a table summarizing all WES analyses and a table for manual' \
                                     + ' QC of WES results.')
    parser.add_argument('--tnp', type=str, help='Path to table of WES tumor/normal pairs.')
    parser.add_argument('--somatic_facets', type=str, nargs="*", help='Optional. Paths to FACETS vcf files.',
                        default=[])
    parser.add_argument('--somatic_maf', type=str, nargs="*", help='Optional. Paths to mutations files.',
                        default=[])
    parser.add_argument('--somatic_maf_agg', type=str, help='Optional. Path to aggregated mutation file.',
                        default=None)
    parser.add_argument('--qc_fastqc_r1', type=str, nargs="*", help='Optional. Paths to FASTQC R1 files.',
                        default=[])
    parser.add_argument('--qc_fastqc_r2', type=str, nargs="*", help='Optional. Paths to FASTQC R2 files.',
                        default=[])
    parser.add_argument('--qc_ncm', type=str, nargs="*", help='Optional. Path to NGSCheckMate files.',
                        default=[])
    parser.add_argument('--qc_flagstat', type=str, nargs="*", help='Optional. Paths to samtools flagstat files.',
                        default=[])
    parser.add_argument('--qc_hsmetrics', type=str, nargs="*", help='Optional. Path to CollectHsMetrics files.',
                        default=[])
    parser.add_argument('--min_cov_t', type=int,
                        help='Optional. Tumor samples with median coverage lower than this value will be flagged.',
                        default=40)
    parser.add_argument('--min_cov_n', type=int,
                        help='Optional. Normal samples with median coverage lower than this value will be flagged.',
                        default=40)
    parser.add_argument('--min_10x_pct_t', type=int,
                        help='Optional. Tumor samples with 10x coverage percent lower than this value will be flagged.',
                        default=60)
    parser.add_argument('--min_10x_pct_n', type=int,
                        help='Optional. Normal samples with 10x coverage percent lower than this value will be flagged.',
                        default=60)
    parser.add_argument('--not_human', action='store_true',
                        help='Optional. Use this option to prevent percentage and million reads formatting in' + \
                        ' human-readable format.')
    parser.add_argument('--output', type=str, help='Path to output table summarizing all WES analyses.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)

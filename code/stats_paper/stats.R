# @created: 25 Nov 24
# @modified: 25 Nov 24
# @authors: Yoann Pradat
#
#    Institut Gustave Roussy
#    Prism Center
#    114 rue Edouard Vaillant, Villejuif, 94800 France
#

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tidyverse))

source("../common/functions/utils.R")

# functions ============================================================================================================

load_mut_status <- function(filepath, samples){
  df_dat_cnt <- load_table(filepath)
  if (is.null(df_dat_cnt)) return(dfs)
  col_row <- "Sample_Id_DNA_T"

  # reformat table
  df_dat_cnt <- df_dat_cnt %>%
    group_by(Sample_Id, Alteration) %>%
    summarize(Count=n()) %>%
    spread(Alteration, Count) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    rename(Sample_Id_DNA_T=Sample_Id) %>%
    mutate(Sample_Id_DNA_T=paste0(Sample_Id_DNA_T, "-DNA")) %>%
    filter(Sample_Id_DNA_T %in% samples)

  # add missing samples if any
  df_dat_cnt <- full_join(df_dat_cnt, tibble(Sample_Id_DNA_T=samples)) %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))

  df_dat_cnt <- df_dat_cnt %>% column_to_rownames(var=col_row)
  df_dat_sts <- df_dat_cnt %>% mutate_if(is.numeric, function(x) as.integer(as.logical(x)))

  df_dat_sts
}


main <- function(args){
  # id of tumor-normal or tumor-NA
  col_p <- "DNA_P"

  # all patients
  df_cln <- read_tsv(args$cln, show_col_types=F, progress=F)
  col_tsb <- "Sample_Id_DNA_T_BAS"
  col_nsb <- "Sample_Id_DNA_N"
  df_cln <- df_cln %>% unite(!!col_p, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

  # load maf file
  df_dna_mut <- load_table(args$dna_mut, header_prefix="##")
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  df_dna_mut <- df_dna_mut %>% unite(!!col_p, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=T)

  # identify ESR1-mutated samples
  gene <- "ESR1"
  df_gene_status <- df_dna_mut %>% filter(Hugo_Symbol=="ESR1") %>% distinct(Hugo_Symbol, DNA_P) %>% mutate(Mutated=1) %>%
    spread(Hugo_Symbol, Mutated)
  df_cln <- left_join(df_cln, df_gene_status)
  df_cln <- df_cln %>% mutate(!!gene := replace_na(!!sym(gene), 0))

  # count total number of mutations
  df_gene_mut_cnt <- df_dna_mut %>% distinct(Hugo_Symbol, DNA_P) %>% group_by(DNA_P) %>% summarize(N_Mut_Genes=n())
  df_cln <- left_join(df_cln, df_gene_mut_cnt)

  # compute TMB
  df_mut_cnt <- df_dna_mut %>% group_by(DNA_P) %>% summarize(N_Mut=n())
  filepath_bed <- "../../data/resources/SureSelect_CR2_padded_10n.bed"
  df_bed <- read_tsv(filepath_bed, col_names=F, show_col_types=F, progress=F)
  target_size <- sum(df_bed$X3-df_bed$X2)/1e6
  df_mut_cnt <- df_mut_cnt %>% rename(TMB="N_Mut") %>% mutate(TMB=TMB/target_size)
  df_cln <- left_join(df_cln, df_mut_cnt)

  # load cna file
  df_dna_cna <- load_table(args$dna_cna)
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  df_dna_cna <- df_dna_cna %>% unite(!!col_p, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=T)

  # count total number of deleted/amplified genes
  df_gene_del_cnt <- df_dna_cna %>% filter(Copy_Number_More %in% c("LOH", "hom_del")) %>%
    group_by(DNA_P) %>% summarize(N_Del_Genes=n())
  df_cln <- left_join(df_cln, df_gene_del_cnt)

  df_gene_amp_cnt <- df_dna_cna %>% filter(Copy_Number_More %in% c("LL_gain", "ML_gain", "HL_gain")) %>%
    group_by(DNA_P) %>% summarize(N_Amp_Genes=n())
  df_cln <- left_join(df_cln, df_gene_amp_cnt)

  # select patients and good-quality baseline WES samples
  df_cln_sub <- df_cln %>% filter(Cancer_Type=="Breast Cancer", QC_DNA_BAS=="PASS")
  df_cln_sub <- df_cln_sub %>% mutate(N_Mut_Genes=replace_na(N_Mut_Genes, 0)) %>%
    mutate(N_Del_Genes=replace_na(N_Del_Genes, 0)) %>%
    mutate(N_Amp_Genes=replace_na(N_Amp_Genes, 0))


  # compute stats
  df_stats <- df_cln_sub %>% group_by(ESR1) %>%
    summarize(Size=n(), Min_TMB=min(TMB), Max_TMB=max(TMB), Median_TMB=median(TMB),
              IQR_TMB=quantile(TMB, 0.75)-quantile(TMB, 0.25),
              Min_N_Mut_Genes=min(N_Mut_Genes),  Max_N_Mut_Genes=max(N_Mut_Genes),
              Median_N_Mut_Genes=median(N_Mut_Genes),
              IQR_N_Mut_Genes=quantile(N_Mut_Genes, 0.75)-quantile(N_Mut_Genes, 0.25),
              Min_N_Del_Genes=min(N_Del_Genes),  Max_N_Del_Genes=max(N_Del_Genes),
              Median_N_Del_Genes=median(N_Del_Genes),
              IQR_N_Del_Genes=quantile(N_Del_Genes, 0.75)-quantile(N_Del_Genes, 0.25),
              Min_N_Amp_Genes=min(N_Amp_Genes),  Max_N_Amp_Genes=max(N_Amp_Genes),
              Median_N_Amp_Genes=median(N_Amp_Genes),
              IQR_N_Amp_Genes=quantile(N_Amp_Genes, 0.75)-quantile(N_Amp_Genes, 0.25))


  write.xlsx(df_stats, args$stats)
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Prepare tables with indices for bootstrapping or cross-validation.')
  parser$add_argument('--cln', type="character", help='Path to data table',
                      default="../../data/cln/curated/cln_icarus_curated_with_qc.tsv")
  parser$add_argument("--bed", type="character", help="Path to BED file",
                      default="../../data/resources/SureSelect_CR2_padded_10n.bed")
  parser$add_argument("--dna_mut", type="character", help="Path to table of all mutations.",
                      default="../../data/wes/somatic_maf/somatic_calls.maf.gz")
  parser$add_argument("--dna_cna", type="character", help="Path to table of all gene cnas.",
                      default="../../data/wes/somatic_cna/somatic_calls.tsv.gz")
  parser$add_argument("--stats", type="character", help="Path to output excel workbook.",
                      default="../../results/stats_paper/stats_breast_baseline.xlsx")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}

# @created: 01 May 22
# @modified: 28 Mar 24
# @authors: Yoann Pradat
#
# Aggregate alterations (copy-number, mutations) annotated with OncoKb and CiVIC annotations.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source("../common/functions/utils.R")

# functions ============================================================================================================

select_samples <- function(df, df_qc_sum){
  col_pair <- "DNA_P"
  df <- df %>% unite(!!col_pair, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_", remove=F)
  mask <- df[[col_pair]] %in% df_qc_sum[[col_pair]]
  n_pair_bef <- length(unique(df[[col_pair]]))
  n_row_bef <- nrow(df)
  df <- df[mask,]
  n_pair_aft <- length(unique(df[[col_pair]]))
  n_row_aft <- nrow(df)
  df <- df %>% select(-all_of(col_pair))
  cat(paste("-INFO: selected data from", paste0(n_pair_aft,"/",n_pair_bef), "QC PASS pairs\n"))
  cat(paste("-INFO: selected", paste0(n_row_aft,"/",n_row_bef), "lines\n"))
  df
}


get_oncokb_columns <- function(){
  c("GENE_IN_ONCOKB", "VARIANT_IN_ONCOKB", "MUTATION_EFFECT", "ONCOGENIC", "TX_CITATIONS",
    "LEVEL_1","LEVEL_2","LEVEL_3A","LEVEL_3B","LEVEL_4","LEVEL_R1","LEVEL_R2","LEVEL_R3", "HIGHEST_LEVEL",
    "LEVEL_Dx1","LEVEL_Dx2","LEVEL_Dx3", "HIGHEST_DX_LEVEL",
    "LEVEL_Px1","LEVEL_Px2","LEVEL_Px3", "HIGHEST_PX_LEVEL")
}


get_civic_columns <- function(){
  c("CIViC_Matching_Disease", "CIViC_Matching_Type", "CIViC_Matching_Gene_Variant",
    "CIViC_Matching_Evidence_Id", "CIViC_Matching_Citation",
    "Predictive:N:A","Predictive:N:B","Predictive:N:C","Predictive:N:D","Predictive:N:E",
    "Predictive:P:A","Predictive:P:B","Predictive:P:C","Predictive:P:D","Predictive:P:E",
    "Diagnostic:N:A","Diagnostic:N:B","Diagnostic:N:C","Diagnostic:N:D","Diagnostic:N:E",
    "Diagnostic:P:A","Diagnostic:P:B","Diagnostic:P:C","Diagnostic:P:D","Diagnostic:P:E",
    "Prognostic:N:A","Prognostic:N:B","Prognostic:N:C","Prognostic:N:D","Prognostic:N:E",
    "Prognostic:P:A","Prognostic:P:B","Prognostic:P:C","Prognostic:P:D","Prognostic:P:E")
}


get_alteration_category_mut <- function(x){
  if (grepl("Ins", x)){
    return("Ins")
  } else if (grepl("Del", x)) {
    return("Del")
  } else {
    return("Mut")
  }
}

get_count_comment_lines <- function(filepath){
  con = file(filepath, "r")
  count <- 0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if (grepl("^#", line)){
      count <- count+1
    } else {
      break
    }
  }

  close(con)
  count
}


process_mut_annotated_line <- function(x){
  x %>%
    gsub(";", "&", .) %>%
    sprintf("?%s", .) %>%
    parse_url() %>%
    extract2("query") %>%
    data.frame(stringsAsFactors=FALSE)
}

process_mut_annotated <- function(df_mut_annotated){
  df_mut_annotated_extra <- bind_rows(sapply(df_mut_annotated$Extra, process_mut_annotated_line))
  df_mut_annotated_extra <- as_tibble(df_mut_annotated_extra)
  df_mut_annotated$Extra <- NULL
  df_mut_annotated <- bind_cols(df_mut_annotated, df_mut_annotated_extra)

  df_mut_annotated %>% separate(`#Uploaded_variation`, c("Uploaded_Feature", "Uploaded_HGVSc"), sep=":", remove=F)
}


aggregate_alterations <- function(df_cna=NULL, df_mut=NULL){
  cat("-aggregating alterations across modalities...")
  cols_civic <- get_civic_columns()
  cols_oncokb <- get_oncokb_columns()
  cols_keep <- c("SAMPLE_ID", "ALTERATION_CATEGORY", cols_oncokb, cols_civic)

  if (!is.null(df_cna)){
    cols_keep_cna <- intersect(colnames(df_cna), c(cols_keep, "HUGO_SYMBOL", "COPY_NUMBER_MORE"))
    df_cna <- df_cna[cols_keep_cna]
  } else {
    df_cna <- data.frame()
  }

  if (!is.null(df_mut)){
    cols_keep_mut <- intersect(colnames(df_mut),
                               c(cols_keep, "HUGO_SYMBOL", "PROTEIN_CHANGE", "PROTEIN_CHANGE_MORE",
                                 "Variant_Classification", "EXON", "t_vaf"))
    df_mut <- df_mut[cols_keep_mut]
  } else {
    df_mut <- data.frame()
  }
  cat("done!\n")

  bind_rows(df_cna, df_mut)
}


make_identifier <- function(df, cols){
  df[cols] <- df[cols] %>% replace(is.na(.), "N/A")
  df <- df %>% unite("Row_Identifier", all_of(cols), remove=F, sep="/")
  df
}


collapse_oncokb_annotations <- function(df_agg){
  cat("-collapsing oncokb annotations...")
  col_cite <- "TX_CITATIONS"
  df_agg <- df_agg %>% rename(Oncokb_Oncogenic=ONCOGENIC, Oncokb_Mutation_Effect=MUTATION_EFFECT)

  #### sensitivity
  levels_sen <- sort(colnames(df_agg)[grepl("^(?=LEVEL_[0-9]).*", colnames(df_agg), perl=T)])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Oncokb_Sen_Level=NA, Oncokb_Sen_Drug=NA, Oncokb_Sen_Citations=NA)
  for (level in levels_sen){
    df_agg_best <- df_agg_best %>%
      mutate(Oncokb_Sen_Level=ifelse(is.na(Oncokb_Sen_Level)&!is.na(.data[[level]]), level, Oncokb_Sen_Level)) %>%
      mutate(Oncokb_Sen_Drug=ifelse(Oncokb_Sen_Level==level, .data[[level]], Oncokb_Sen_Drug)) %>%
      mutate(Oncokb_Sen_Citations=ifelse(Oncokb_Sen_Level==level, .data[[col_cite]], Oncokb_Sen_Citations))
  }
  df_agg_best <- df_agg_best %>% mutate(Oncokb_Sen_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Oncokb_Sen_Level=NA, Oncokb_Sen_Drug=NA, Oncokb_Sen_Citations=NA,
                                  Oncokb_Oncogenic=NA, Oncokb_Mutation_Effect=NA)
  dfs_agg_all <- list()
  for (level in levels_sen){
    df_agg_level <- df_agg_all %>%
      mutate(Oncokb_Sen_Level=level) %>%
      mutate(Oncokb_Sen_Drug=.data[[level]]) %>%
      mutate(Oncokb_Sen_Citations=.data[[col_cite]]) %>%
      filter(!is.na(Oncokb_Sen_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Oncokb_Sen_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Oncokb_Sen_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Oncokb_Sen_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  sen_level_simple <- list(LEVEL_1="Tier1", LEVEL_2="Tier1", LEVEL_3A="Tier2", LEVEL_3B="Tier3", LEVEL_4="Tier3")
  if (!all(is.na(df_agg$Oncokb_Sen_Level))){
    df_agg$Oncokb_Sen_Level_Simple <- recode(df_agg$Oncokb_Sen_Level, !!!sen_level_simple)
  } else {
    df_agg$Oncokb_Sen_Level_Simple <- NA
  }


  #### resistance
  levels_res <- sort(colnames(df_agg)[grepl("^LEVEL_R", colnames(df_agg))])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Oncokb_Res_Level=NA, Oncokb_Res_Drug=NA, Oncokb_Res_Citations=NA)
  for (level in levels_res){
    df_agg_best <- df_agg_best %>%
      mutate(Oncokb_Res_Level=ifelse(is.na(Oncokb_Res_Level)&!is.na(.data[[level]]), level, Oncokb_Res_Level)) %>%
      mutate(Oncokb_Res_Drug=ifelse(Oncokb_Res_Level==level, .data[[level]], Oncokb_Res_Drug)) %>%
      mutate(Oncokb_Res_Citations=ifelse(Oncokb_Res_Level==level, .data[[col_cite]], Oncokb_Res_Citations))
  }
  df_agg_best <- df_agg_best %>% mutate(Oncokb_Res_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Oncokb_Res_Level=NA, Oncokb_Res_Drug=NA, Oncokb_Res_Citations=NA)
  dfs_agg_all <- list()
  for (level in levels_res){
    df_agg_level <- df_agg_all %>%
      mutate(Oncokb_Res_Level=level) %>%
      mutate(Oncokb_Res_Drug=.data[[level]]) %>%
      mutate(Oncokb_Res_Citations=.data[[col_cite]]) %>%
      filter(!is.na(Oncokb_Res_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Oncokb_Res_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Oncokb_Res_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Oncokb_Res_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  res_level_simple <- list(LEVEL_R1="Tier1", LEVEL_R2="Tier2", LEVEL_R3="Tier3")
  if (!all(is.na(df_agg$Oncokb_Res_Level))){
    df_agg$Oncokb_Res_Level_Simple <- recode(df_agg$Oncokb_Res_Level, !!!res_level_simple)
  } else {
    df_agg$Oncokb_Res_Level_Simple <- NA
  }

  cat("done!\n")
  cols <- c("Oncokb_Annotated", "Oncokb_Oncogenic", "Oncokb_Mutation_Effect",
            "Oncokb_Sen_Level", "Oncokb_Sen_Level_Simple", "Oncokb_Sen_Drug", "Oncokb_Sen_Drug_Best_Level",
            "Oncokb_Sen_Citations", "Oncokb_Res_Level", "Oncokb_Res_Level_Simple",
            "Oncokb_Res_Drug", "Oncokb_Res_Drug_Best_Level", "Oncokb_Res_Citations", "Row_Id")
  df_agg %>% select(all_of(cols))
}


collapse_civic_annotations <- function(df_agg){
  cat("-collapsing civic annotations...")
  col_cite <- "CIViC_Matching_Citation"
  col_evid <- "CIViC_Matching_Evidence_Id"

  #### sensitivity
  levels_sen <- sort(colnames(df_agg)[grepl("^Predictive:P", colnames(df_agg), perl=T)])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Civic_Sen_Level=NA, Civic_Sen_Drug=NA, Civic_Sen_Citations=NA)
  for (level in levels_sen){
    df_agg_best <- df_agg_best %>%
      mutate(Civic_Sen_Level=ifelse(is.na(Civic_Sen_Level)&!is.na(.data[[level]]), level, Civic_Sen_Level)) %>%
      mutate(Civic_Sen_Drug=ifelse(Civic_Sen_Level==level, .data[[level]], Civic_Sen_Drug)) %>%
      mutate(Civic_Sen_Citations=ifelse(Civic_Sen_Level==level, .data[[col_cite]], Civic_Sen_Citations)) %>%
      mutate(Civic_Sen_Evidence_Id=ifelse(Civic_Sen_Level==level, .data[[col_evid]], Civic_Sen_Evidence_Id))
  }
  df_agg_best <- df_agg_best %>% mutate(Civic_Sen_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Civic_Sen_Level=NA, Civic_Sen_Drug=NA, Civic_Sen_Citations=NA,
                                  Civic_Sen_Evidence_Id=NA)
  dfs_agg_all <- list()
  for (level in levels_sen){
    df_agg_level <- df_agg_all %>%
      mutate(Civic_Sen_Level=level) %>%
      mutate(Civic_Sen_Drug=.data[[level]]) %>%
      mutate(Civic_Sen_Citations=.data[[col_cite]]) %>%
      mutate(Civic_Sen_Evidence_Id=.data[[col_evid]]) %>%
      filter(!is.na(Civic_Sen_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Civic_Sen_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Civic_Sen_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Civic_Sen_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  sen_level_simple <- list(`Predictive:P:A`="Tier1", `Predictive:P:B`="Tier2", `Predictive:P:C`="Tier3",
                           `Predictive:P:D`="Tier3", `Predictive:P:E`="Tier3")
  if (!all(is.na(df_agg$Civic_Sen_Level))){
    df_agg$Civic_Sen_Level_Simple <- recode(df_agg$Civic_Sen_Level, !!!sen_level_simple)
  } else {
    df_agg$Civic_Sen_Level_Simple <- NA
  }


  #### resistance
  levels_res <- sort(colnames(df_agg)[grepl("^Predictive:N", colnames(df_agg))])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Civic_Res_Level=NA, Civic_Res_Drug=NA, Civic_Res_Citations=NA,
                                   Civic_Res_Evidence_Id=NA)
  for (level in levels_res){
    df_agg_best <- df_agg_best %>%
      mutate(Civic_Res_Level=ifelse(is.na(Civic_Res_Level)&!is.na(.data[[level]]), level, Civic_Res_Level)) %>%
      mutate(Civic_Res_Drug=ifelse(Civic_Res_Level==level, .data[[level]], Civic_Res_Drug)) %>%
      mutate(Civic_Res_Citations=ifelse(Civic_Res_Level==level, .data[[col_cite]], Civic_Res_Citations)) %>%
      mutate(Civic_Res_Evidence_Id=ifelse(Civic_Res_Level==level, .data[[col_evid]], Civic_Res_Evidence_Id))
  }
  df_agg_best <- df_agg_best %>% mutate(Civic_Res_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Civic_Res_Level=NA, Civic_Res_Drug=NA, Civic_Res_Drug=NA, Civic_Res_Citations=NA,
                                  Civic_Res_Evidence_Id=NA)
  dfs_agg_all <- list()
  for (level in levels_res){
    df_agg_level <- df_agg_all %>%
      mutate(Civic_Res_Level=level) %>%
      mutate(Civic_Res_Drug=.data[[level]]) %>%
      mutate(Civic_Res_Citations=.data[[col_cite]]) %>%
      mutate(Civic_Res_Evidence_Id=.data[[col_evid]]) %>%
      filter(!is.na(Civic_Res_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Civic_Res_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Civic_Res_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Civic_Res_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  res_level_simple <- list(`Predictive:N:A`="Tier1", `Predictive:N:B`="Tier2", `Predictive:N:C`="Tier3",
                           `Predictive:N:D`="Tier3", `Predictive:N:E`="Tier3")
  if (!all(is.na(df_agg$Civic_Res_Level))){
    df_agg$Civic_Res_Level_Simple <- recode(df_agg$Civic_Res_Level, !!!res_level_simple)
  } else {
    df_agg$Civic_Res_Level_Simple <- NA
  }


  cat("done!\n")
  cols <- c("Civic_Annotated",
            "Civic_Sen_Level", "Civic_Sen_Level_Simple", "Civic_Sen_Drug", "Civic_Sen_Drug_Best_Level",
            "Civic_Sen_Citations", "Civic_Sen_Evidence_Id", "Civic_Res_Level", "Civic_Res_Level_Simple",
            "Civic_Res_Drug", "Civic_Res_Drug_Best_Level", "Civic_Res_Citations", "Civic_Res_Evidence_Id", "Row_Id")
  df_agg %>% select(all_of(cols))
}


get_alteration <- function(df_agg){
  regex <- "(?<=\\/)[0-9]*"
  cols <- c("Hugo_Symbol", "Alteration_Category", "Alteration", "Alteration_Detail")
  df_agg %>%
    mutate(Alteration=NA) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Deletion", "Del", Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Amplification", "Amp", Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Fusion", "Fus", Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Mut",gsub("^p.","",PROTEIN_CHANGE),Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Ins",
                             paste("Exon", gsub(regex,"",EXON,perl=T), "Ins"), Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Del",
                             paste("Exon", gsub(regex,"",EXON,perl=T), "Del"), Alteration)) %>%
    mutate(Alteration_Detail=ifelse(ALTERATION_CATEGORY=="Mut",PROTEIN_CHANGE,NA)) %>%
    mutate(Alteration_Detail=ifelse(ALTERATION_CATEGORY=="Ins",PROTEIN_CHANGE,Alteration_Detail)) %>%
    mutate(Alteration_Detail=ifelse(ALTERATION_CATEGORY=="Del",PROTEIN_CHANGE,Alteration_Detail)) %>%
    mutate(Alteration_Detail=ifelse(!is.na(COPY_NUMBER_MORE), COPY_NUMBER_MORE, Alteration_Detail)) %>%
    mutate(Alteration=ifelse(!is.na(PROTEIN_CHANGE_MORE), PROTEIN_CHANGE_MORE, Alteration)) %>%
    mutate(Alteration_Detail=ifelse(!is.na(PROTEIN_CHANGE_MORE), PROTEIN_CHANGE, Alteration_Detail)) %>%
    mutate(Alteration=ifelse(!is.na(Alteration), gsub("\\/","",Alteration,perl=T), NA)) %>%
    mutate(Alteration=ifelse(is.na(Alteration), Variant_Classification, Alteration)) %>%
    mutate(HUGO_SYMBOL=ifelse(ALTERATION_CATEGORY=="Fusion", FUSION, HUGO_SYMBOL)) %>%
    rename(Alteration_Category=ALTERATION_CATEGORY, Hugo_Symbol=HUGO_SYMBOL) %>%
    select(all_of(cols))
}


harmonize_drugs <- function(x, df_drug){
  x <- x[!is.na(x)]
  if (length(x)==0){
    return(NA)
  } else {
    vals <- unlist(strsplit(x, ",|;|\\+|\\|"))
    vals <- sapply(vals, function(v) toupper(trimws(v)))
    vals <- toupper(vals)
    if (!all(vals %in% df_drug$Drug)){
      vals_mis <- setdiff(vals, df_drug$Drug)
      warning(paste("-the following drugs are not in the drug table:", paste(vals_mis, collapse=",")))
    } else {
      vals_mis <- c()
    }

    x <- toupper(x)
    df_recode <- df_drug %>% filter(Drug %in% vals) %>% select(Drug, DCI) %>% distinct()
    df_recode <- bind_rows(df_recode, tibble(Drug=vals_mis, DCI=vals_mis))
    for (i in 1:nrow(df_recode)){
      x <- gsub(df_recode[i,"Drug"], df_recode[i,"DCI"], x)
    }

    vals <- unlist(strsplit(x, ",|;"))
    vals <- sapply(vals, function(v) toupper(trimws(v)))
    vals <- sort(unique(vals))
    return(paste0(vals, collapse="|"))
  }
}


union_drugs <- function(vals){
  vals <- vals[!is.na(vals)]
  if (length(vals)==0){
    return(NA)
  } else {
    vals <- paste0(vals, collapse=",")
    vals <- unlist(strsplit(vals, ",|;"))
    vals <- sapply(vals, function(v) toupper(trimws(v)))
    vals <- sort(unique(vals))
    return(paste0(vals, collapse="|"))
  }
}


save_table <- function(table, output){
  if (grepl(".gz$", output)){
    write.table(table, file=gsub(".gz$", "", output), sep="\t", quote=F, row.names=F, na="")
    system(paste("gzip", gsub(".gz$", "", output)))
  } else {
    write.table(table, file=output, sep="\t", quote=F, row.names=F, na="")
  }
  cat(paste("-file saved at", output, "\n"))
}


main <- function(args){
  # load tables
  df_qc_sum <- load_table(args$qc_summary)
  df_ids <- load_table(args$ids)
  df_cln <- load_table(args$cln) %>% rename(Tumor_Type=Cancer_Type)
  df_cna <- load_table(args$cna, guess_max=1e5)
  df_cna_pass <- load_table(args$cna_pass, guess_max=1e5)
  df_mut <- load_table(args$mut, guess_max=1e5)

  # select samples
  df_qc_sum <- df_qc_sum %>% filter(QC_Final_Decision=="PASS")
  df_cna <- select_samples(df_cna, df_qc_sum)
  df_cna_pass <- select_samples(df_cna_pass, df_qc_sum)
  df_mut <- select_samples(df_mut, df_qc_sum)

  # process cnas
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_psb <- "Pair_Sample_Barcode"
  df_cna <- df_cna %>% mutate(ALTERATION_CATEGORY=Alteration, SAMPLE_ID=Tumor_Sample_Barcode, HUGO_SYMBOL=Hugo_Symbol)
  df_cna_pass <- df_cna_pass %>% mutate(COPY_NUMBER_MORE=paste(Copy_Number_More, "-", `TCN_EM:LCN_EM`))
  df_cna <- df_cna %>% unite(!!col_psb, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
  df_cna_pass <- df_cna_pass %>% unite(!!col_psb, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
  df_cna <- left_join(df_cna, df_cna_pass %>% select(all_of(col_tsb), Hugo_Symbol, COPY_NUMBER_MORE),
                      by=c(col_tsb, "Hugo_Symbol"))
  df_cna <- df_cna %>% select(-all_of(col_psb))

  # process mutations:
  # - where TP53, replace PROTEIN_CHANGE by DNA Binding Domain where applicable
  df_mut <- df_mut %>% rowwise %>% mutate(ALTERATION_CATEGORY=get_alteration_category_mut(Variant_Classification)) %>%
    rename(SAMPLE_ID=Tumor_Sample_Barcode, HUGO_SYMBOL=Hugo_Symbol, PROTEIN_CHANGE=HGVSp_Short) %>%
    mutate(PROTEIN_CHANGE=ifelse(is.na(PROTEIN_CHANGE), HGVSc, PROTEIN_CHANGE))
  df_mut[["PROTEIN_CHANGE"]] <- gsub("%3D", "=", df_mut[["PROTEIN_CHANGE"]])
  mask_tp53 <- df_mut[["HUGO_SYMBOL"]]=="TP53"
  mask_dnabd <- (df_mut$Start_Position >= 7577149) & (df_mut$End_Position <- 7578443)
  df_mut[mask_tp53 & mask_dnabd, "PROTEIN_CHANGE_MORE"] <- "DNA Binding Domain"

  mask_splice <- grepl("splice", df_mut[["PROTEIN_CHANGE"]])
  mask_splice <- mask_splice | grepl("Splice|Silent", df_mut[["Variant_Classification"]])
  df_mut[mask_splice, "PROTEIN_CHANGE_MORE"] <- "Splice_Site"
  df_mut[mask_splice, "Variant_Classification"] <- "Splice_Site"
  df_mut["t_vaf"] = df_mut["t_alt_count"]/df_mut["t_depth"]

  # aggregate
  df_agg <- aggregate_alterations(df_cna=df_cna, df_mut=df_mut)
  df_agg <- df_agg %>% distinct()
  df_agg_alteration <- get_alteration(df_agg)

  # add indicators of oncokb and civic annotation
  mask_okb <- !is.na(df_agg["ONCOGENIC"])
  mask_civ <- !is.na(df_agg["CIViC_Matching_Gene_Variant"])
  df_agg[mask_okb, "Oncokb_Annotated"] <- "Yes"
  df_agg[mask_civ, "Civic_Annotated"] <- "Yes"

  # add row id because collapse of oncokb and civic annotations may produce multiple
  # rows per row id in case multiple levels of annotation coexist for the same alteration
  df_agg <- df_agg %>% mutate(Row_Id=1:nrow(df_agg))

  # collapse annotations
  df_agg_oncokb <- collapse_oncokb_annotations(df_agg)
  df_agg_civic <- collapse_civic_annotations(df_agg)

  # add Subject_Id and Tumor_Type columns
  df_agg <- df_agg %>% rename(Sample_Id=SAMPLE_ID)
  cols_ids <- c("Sample_Id", "Biopsy_Id", "Subject_Id", "Biopsy_Visit")
  df_agg <- left_join(df_agg, df_ids[,cols_ids] %>% distinct(), by="Sample_Id")
  cols_cln <- c("Subject_Id", "Tumor_Type", "MSKCC_Oncotree", "Civic_Disease")
  df_agg <- left_join(df_agg, df_cln[,cols_cln] %>% distinct(), by="Subject_Id")

  # column bind and joins
  cols_cln <- setdiff(cols_cln, cols_ids)
  cols_agg <- c("Variant_Classification", "t_vaf")
  df_fin <- bind_cols(df_agg[c(cols_ids, cols_cln, cols_agg, "Row_Id")], df_agg_alteration)
  cols_fin <- colnames(df_fin)
  df_fin <- make_identifier(df_fin, setdiff(colnames(df_fin), c("Row_Id", "t_vaf")))

  df_agg_oncokb_best <- df_agg_oncokb %>% filter(Oncokb_Sen_Drug_Best_Level==1, Oncokb_Res_Drug_Best_Level==1)
  df_agg_oncokb_oth <- df_agg_oncokb %>% filter(Oncokb_Sen_Drug_Best_Level!=1 | Oncokb_Res_Drug_Best_Level!=1)
  df_agg_civic_best <- df_agg_civic %>% filter(Civic_Sen_Drug_Best_Level==1, Civic_Res_Drug_Best_Level==1)
  df_agg_civic_oth <- df_agg_civic %>% filter(Civic_Sen_Drug_Best_Level!=1 | Civic_Res_Drug_Best_Level!=1)

  df_fin_best <- full_join(df_fin, df_agg_oncokb_best, by="Row_Id")
  df_fin_best <- full_join(df_fin_best, df_agg_civic_best, by="Row_Id")

  df_fin_oth_oncokb <- right_join(df_fin, df_agg_oncokb_oth, by="Row_Id")
  df_fin_oth_civic <- right_join(df_fin, df_agg_civic_oth, by="Row_Id")
  df_fin_oth <- full_join(df_fin_oth_oncokb, df_fin_oth_civic, relationship="many-to-many")
  df_fin_all <- bind_rows(df_fin_best, df_fin_oth) %>% arrange(Row_Identifier)

  df_fin_all$Row_Id <- NULL
  df_fin_all <- df_fin_all %>% distinct()

  # for CNA, set Variant_Classification to Deletion or Amplification
  mask_null <- is.na(df_fin_all$Variant_Classification)
  df_fin_all[mask_null, "Variant_Classification"] <- df_fin_all[mask_null, "Alteration_Category"]

  # consensus level by taking the best of oncokb and civic
  cols <- c("Oncokb_Sen_Level_Simple", "Civic_Sen_Level_Simple")
  df_fin_all$Sen_Level_Simple <- pmin(df_fin_all[[cols[1]]], df_fin_all[[cols[2]]], na.rm=T)
  cols <- c("Oncokb_Res_Level_Simple", "Civic_Res_Level_Simple")
  df_fin_all$Res_Level_Simple <- pmin(df_fin_all[[cols[1]]], df_fin_all[[cols[2]]], na.rm=T)

  # harmonize drug names
  df_drug <- load_table(args$drug)
  cols_drug <- c("Oncokb_Sen_Drug", "Oncokb_Res_Drug", "Civic_Sen_Drug", "Civic_Res_Drug")
  for (col_drug in cols_drug){
    df_uni <- tibble(Old=unique(df_fin_all[[col_drug]]))
    df_uni$New <- unlist(lapply(df_uni$Old, function(x) harmonize_drugs(x, df_drug)))
    df_fin_all <- left_join(df_fin_all, df_uni, by=setNames(c("Old"), col_drug))
    df_fin_all <- df_fin_all %>% select(-all_of(col_drug)) %>% rename(!!col_drug:=New)
  }

  # select only annotations with best level
  df_fin_best <- df_fin_all %>% filter_at(vars(ends_with("_Best_Level")), all_vars(. == 1))

  # save tables
  save_table(df_fin_all, args$output_all)
  save_table(df_fin_best, args$output_best)
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Aggregate alterations.')
  parser$add_argument("--qc_summary", type="character", help="Path to table of QC summary.",
                      default="../../data/wes/wes_analyses_summary_b1_b2_b3_b4.xlsx")
  parser$add_argument("--ids", type="character", help="Path to input curated ids table.",
                      default="../../data/cln/curated/ids_icarus_curated.tsv")
  parser$add_argument("--cln", type="character", help="Path to input curated clinical table.",
                      default="../../data/cln/curated/cln_icarus_curated.tsv")
  parser$add_argument("--cna", type="character", help="Path to input annotated CNAs table.",
                      default="../../data/wes/somatic_cna/somatic_calls_union_ann.tsv.gz")
  parser$add_argument("--cna_pass", type="character", help="Path to input annotated CNAs pass table.",
                      default="../../data/wes/somatic_cna/somatic_calls.tsv.gz")
  parser$add_argument("--mut", type="character", help="Path to input annotated mutations table.",
                      default="../../data/wes/somatic_maf/somatic_calls_union_ann.maf.gz")
  parser$add_argument('--drug', type="character", help='Path to table of drugs.',
                      default="../../../../completed/MetaPRISM/data/resources/drug_tables/Table_Drugs_v8.xlsx")
  parser$add_argument("--output_best", type="character",
                      help="Path to output aggregated table with only the best annotation.",
                      default="../../results/wes_analysis/alterations/aggregated_alterations_best.tsv")
  parser$add_argument("--output_all", type="character",
                      help="Path to output aggregated table with all annotations.",
                      default="../../results/wes_analysis/alterations/aggregated_alterations_all.tsv")
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}

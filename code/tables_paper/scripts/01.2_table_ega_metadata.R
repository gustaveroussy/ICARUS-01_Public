# @created: 16 Dec 22
# @modified: 20 Nov 24
# @authors: Yoann Pradat
#
# Prepare the metadata table for EGA upload.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(openxlsx))

# functions ============================================================================================================

firstdown <- function(x) {
  substr(x, 1, 1) <- tolower(substr(x, 1, 1))
  x
}


# load
load_samples_data <- function(cln_table, sam_wes_table, sam_rna_table, mode="breast"){
  cat("-loading clinical and bio table: ... ")
  cln <- read.delim(cln_table, sep="\t")
  sam_wes <- read.delim(sam_wes_table, sep="\t")
  sam_rna <- read.delim(sam_rna_table, sep="\t")
  cat("ok!\n")

  # mode-specific variables
  if (mode=="breast"){
    cancer_type <- "Breast Cancer"
    design_name_wes <- "ICARUS-BREAST01_WES_NovaSeq_6000_SureSelect_CRE2"
    design_name_rna <- "ICARUS-BREAST01_RNAseq_NovaSeq_6000_SureSelect_Strand_Specific"
  } else if (mode=="lung") {
    cancer_type <- "Lung Cancer"
    design_name_wes <- "ICARUS-LUNG01_WES_NovaSeq_6000_SureSelect_CRE2"
    design_name_rna <- "ICARUS-LUNG01_RNAseq_NovaSeq_6000_SureSelect_Strand_Specific"
  }

  # wes samples
  wes_samples_blood <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_DNA_BAS=="PASS") %>%
    filter(Sample_Id_DNA_N!="") %>%
    pull(var=Sample_Id_DNA_N)
  wes_bio_sites_blood <- rep("Blood", length(wes_samples_blood))

  wes_samples_bas <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_DNA_BAS=="PASS") %>%
    pull(var=Sample_Id_DNA_T_BAS)
  wes_bio_sites_bas <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_DNA_BAS=="PASS") %>%
    pull(var=Baseline_Biopsy_Site_Precise)

  wes_samples <- c(wes_samples_blood, wes_samples_bas)
  wes_bio_sites <- c(wes_bio_sites_blood, wes_bio_sites_bas)

  # rna samples
  rna_samples_bas <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_RNA_BAS=="PASS", QC_RNA_ONT=="PASS") %>%
    pull(var=Sample_Id_RNA_T_BAS)
  rna_bio_sites_bas <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_RNA_BAS=="PASS", QC_RNA_ONT=="PASS") %>%
    pull(var=Baseline_Biopsy_Site_Precise)

  rna_samples_ont <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_RNA_BAS=="PASS", QC_RNA_ONT=="PASS") %>%
    pull(var=Sample_Id_RNA_T_ONT)
  rna_bio_sites_ont <- cln %>%
    filter(Cancer_Type==cancer_type) %>%
    filter(QC_RNA_BAS=="PASS", QC_RNA_ONT=="PASS") %>%
    pull(var=On_treatment_Biopsy_Site)

  rna_samples <- c(rna_samples_bas, rna_samples_ont)
  rna_bio_sites <- c(rna_bio_sites_bas, rna_bio_sites_ont)

  # biopsy table - one line per sample
  bio_wes <- tibble(Sample_Id=wes_samples, Biopsy_Site=wes_bio_sites) %>%
    mutate(Subject_Id=substr(Sample_Id, 1, 8)) %>%
    mutate(Capture_Kit="SureSelect XT Human All Exon CRE-V2") %>%
    mutate(Library_Source="GENOMIC") %>%
    mutate(Library_Selection="Hybrid Selection") %>%
    mutate(Library_Strategy="WXS") %>%
    mutate(Timepoint=ifelse(grepl("-BAS", Sample_Id), "BAS", ifelse(grepl("-N", Sample_Id), "N", "EOT"))) %>%
    mutate(Sample_Type=ifelse(grepl("-N-DNA$", Sample_Id), "DNA_N", "DNA_T")) %>%
    mutate(Design_Name=design_name_wes)
  bio_wes <- left_join(bio_wes, sam_wes %>% select(Sample_Id, FASTQ_1, FASTQ_2))

  bio_rna <- tibble(Sample_Id=rna_samples, Biopsy_Site=rna_bio_sites) %>%
    mutate(Subject_Id=substr(Sample_Id, 1, 8)) %>%
    mutate(Capture_Kit="SureSelect Strand-Specific RNA Library Preparation kit") %>%
    mutate(Library_Source="TRANSCRIPTOMIC") %>%
    mutate(Library_Selection="PolyA") %>%
    mutate(Library_Strategy="RNA-Seq") %>%
    mutate(Timepoint=ifelse(grepl("-BAS", Sample_Id), "BAS", "ONT")) %>%
    mutate(Sample_Type="RNA_T") %>%
    mutate(Design_Name=design_name_rna)
  bio_rna <- left_join(bio_rna, sam_rna %>% select(Sample_Id, FASTQ_1, FASTQ_2))

  # join wes & rna
  bio <- bind_rows(bio_wes, bio_rna) %>%
    mutate(Platform="Illumina NovaSeq 6000")

  # join samples & clinics
  if (mode=="breast"){
    cols_cln <- c("Subject_Id", "Gender", "Diagnosis_HER2_Status", "Diagnosis_HR_Status", "Diagnosis_Histology")
    bio <- bio %>% left_join(cln[, cols_cln], by="Subject_Id")
    bio <- bio %>%
      mutate(Diagnosis_HR_Status = paste("HR", tolower(Diagnosis_HR_Status))) %>%
      mutate(Diagnosis_HER2_Status = ifelse(Diagnosis_HER2_Status == "", "HER2 unknown", Diagnosis_HER2_Status))

    # add phenotype
    bio$Tumor_Type <- "Breast cancer"
    bio <- bio %>% unite(IHC, Diagnosis_HER2_Status, Diagnosis_HR_Status, sep=" ")
    bio <- bio %>% unite(Phenotype, Tumor_Type, Diagnosis_Histology, IHC, sep=", ")
  } else if (mode=="lung") {
    cols_cln <- c("Subject_Id", "Gender", "Diagnosis_Histology")


    # add phenotype
    bio$Tumor_Type <- "Lung cancer"
    bio <- bio %>% left_join(cln[, cols_cln], by="Subject_Id")
    bio <- bio %>% unite(Phenotype, Tumor_Type, Diagnosis_Histology, sep=", ")
  }

  # total files
  n_row_cur <- bio %>% nrow()
  cat(paste("-total number of samples:", n_row_cur), "\n")

  bio
}


print_statistics <- function(bio){
  stats <- list()

  # distinct patients
  nb <- bio %>%
    select (Subject_Id) %>%
    unique() %>%
    nrow()
  stats["total_subjects"] <- nb
  cat("-# of distinct subjects: ... ", nb, "\n")

  # distinct tumor DNA samples in dataset
  nb <- bio %>% filter(Sample_Type == "DNA_T") %>%
    select (Sample_Id) %>%
    unique() %>%
    nrow()
  stats["tumor_dna_samples"] <- nb
  cat("-# of distinct tumor DNA samples: ... ", nb, "\n")

  # distinct blood samples in dataset
  nb <- bio %>% filter(Sample_Type == "DNA_N") %>%
    select (Sample_Id) %>%
    unique() %>%
    nrow()
  stats["normal_dna_samples"] <- nb
  cat("-# of distinct blood DNA samples: ... ", nb, "\n")

  # tumor DNA samples at BAS
  nb <- bio %>% filter(Sample_Type == "DNA_T", Timepoint == "BAS") %>%
    select (Sample_Id) %>%
    unique() %>%
    nrow()
  stats["bas_tumor_dna_samples"] <- nb
  cat("-# of distinct BAS tumor DNA samples: ... ", nb, "\n")

  # distinct tumor RNA samples in dataset
  nb <- bio %>% filter(Sample_Type == "RNA_T") %>%
    select (Sample_Id) %>%
    unique() %>%
    nrow()
  stats["tumor_rna_samples"] <- nb
  cat("-# of distinct tumor RNA samples: ... ", nb, "\n")

  # tumor RNA samples at BAS
  nb <- bio %>% filter(Sample_Type == "RNA_T", Timepoint == "BAS") %>%
    select (Sample_Id) %>%
    unique() %>%
    nrow()
  stats["bas_tumor_rna_samples"] <- nb
  cat("-# of distinct BAS tumor RNA samples: ... ", nb, "\n")

  # tumor RNA samples at ONT
  nb <- bio %>% filter(Sample_Type == "RNA_T", Timepoint == "ONT") %>%
    select (Sample_Id) %>%
    unique() %>%
    nrow()
  stats["ont_tumor_rna_samples"] <- nb
  cat("-# of distinct ONT tumor RNA samples: ... ", nb, "\n")

  stats
}


add_anonymous_identifiers <- function(bio, ids_anonym){
  bio <- left_join(bio, ids_anonym)
  bio <- bio %>%
    rowwise() %>%
    mutate(bioSampleId = gsub(Subject_Id, Subject_Id_Anonymous, Sample_Id)) %>%
    mutate(subjectId = Subject_Id_Anonymous) %>%
    ungroup()

  # some sample ids hold extra info about tube number to discriminate between double biopsies
  # reomve this extra number
  bio <- bio %>%
    rowwise() %>%
    mutate(
      bioSampleId = if_else(
        grepl(paste0("^", subjectId, "-[0-9]+"), bioSampleId),
        gsub(paste0("(", subjectId, ")-[0-9]+"), "\\1", bioSampleId),
        bioSampleId
      )
    )

  stopifnot(sum(is.na(bio$subjectId))==0)
  stopifnot(sum(is.na(bio$bioSampleId))==0)
  cat("-added anonymous identifiers\n")

  bio
}


make_ega_metadata <- function(bio, stats, mode){
  # output metadata
  # title	alias	description	subjectId	bioSampleId	caseOrControl	gender	organismPart	cellLine	region	phenotype

  if (mode=="breast"){
    EGA_title <- paste("Efficacy, safety and biomarker analysis of ICARUS-BREAST01: a phase 2 Study of Patritumab",
                       "Deruxtecan in patients with HR+/HER2- advanced breast cancer")
    EGA_alias <- "ICARUS-BREAST01"
    EGA_description <- paste("A phase II clinical trial (NCT04965766) of patritumab deruxtecan in 99 breast cancer",
                             "patients. Whole-exome sequencing (WES) is available for", stats$tumor_dna_samples,
                             "tumor samples and", stats$normal_dna_samples, "blood samples collected at entry into the",
                             "trial. RNA-squencing is available for", stats$tumor_rna_samples, "samples comprising",
                             stats$bas_tumor_rna_samples, "samples", "collected at entry into the trial and",
                             stats$ont_tumor_rna_samples, "samples collected during treatment (cycle 1 day 3,",
                             "cycle 1 day 19, or cycle 2 day 3) from", stats$bas_tumor_rna_samples, "patients.")
  } else if (mode=="lung") {
    EGA_title <- paste("Efficacy, safety and biomarker analysis of ICARUS-LUNG01: a phase 2 Study of Datopotomab",
                       "Deruxtecan (Dato-DXd) in advanced Non-Small Cell Lung Cancer (NSCLC) patients")
    EGA_alias <- "ICARUS-LUNG01"
    EGA_description <- paste("A phase II clinical trial (NCT04940325) of datopotomab deruxtecan in 100 lung cancer",
                             "patients. Whole-exome sequencing (WES) is available for", stats$tumor_dna_samples,
                             "tumor samples and", stats$normal_dna_samples, "blood samples collected at entry into the",
                             "trial. RNA-squencing is available for", stats$tumor_rna_samples, "samples comprising",
                             stats$bas_tumor_rna_samples, "samples", "collected at entry into the trial and",
                             stats$ont_tumor_rna_samples, "samples collected during treatment (cycle 1 day 3,",
                             "or cycle 2 day 3) from", stats$bas_tumor_rna_samples, "patients.")
  }

  EGA_metadata <- bio %>%
    mutate(title=EGA_title) %>%
    mutate(alias=EGA_alias) %>%
    mutate(description=EGA_description) %>%
    mutate(caseOrControl=tolower(ifelse(Sample_Type=="DNA_N", "Control", "Case"))) %>%
    mutate(gender=tolower(Gender)) %>%
    mutate(organismPart=Biopsy_Site) %>%
    mutate(cellLine="") %>%
    mutate(region="") %>%
    mutate(phenotype=Phenotype) %>%
    mutate(designName=Design_Name) %>%
    mutate(Instrument_Model=Platform) %>%
    mutate(Library_Name=Capture_Kit) %>%
    mutate(Library_Layout="PAIRED") %>%
    mutate(Library_Source=Library_Source) %>%
    mutate(Library_Selection=Library_Selection) %>%
    mutate(Library_Strategy=Library_Strategy) %>%
    mutate(Library_Construction_Protocol="") %>%
    rename(`First Fastq File`=FASTQ_1, `Second Fastq File`=FASTQ_2) %>%
    mutate(`First Unencrypted checksum`="", `Second Unencrypted checksum`="") %>%
    mutate(`First Checksum`="", `Second Checksum`="") %>%
    select(title, alias, description, subjectId, bioSampleId, caseOrControl, gender,
           organismPart, cellLine, region, phenotype, designName, Design_Name,
           Instrument_Model, Library_Name, Library_Layout, Library_Source, Library_Selection,
           Library_Strategy, Library_Construction_Protocol, `First Fastq File`, `First Checksum`,
           `First Unencrypted checksum`, `Second Fastq File`,`Second Checksum`,
           `Second Unencrypted checksum`, Sample_Id)

  # list of sheets
  dfs_ega <- list()

  # sheet Title_Description
  dfs_ega[["Title_Description"]] <- tibble(Submission_Title=EGA_alias, Submission_Description=EGA_alias)

  # sheet Study
  dfs_ega[["Study"]] <- tibble(Study_Descriptive_Title=EGA_title,
                               Study_Type="Cancer Genomics, Exome Sequencing, RNA Sequencing",
                               Short_Name=EGA_alias,
                               Abstract=EGA_description)

  # sheet Samples
  dfs_ega[["Samples"]] <- EGA_metadata %>% select(designName, title, alias, description, subjectId, bioSampleId,
    caseOrControl, gender, organismPart, cellLine, region, phenotype) %>% arrange(designName)

  # sheet Experiments
  dfs_ega[["Experiments"]] <- EGA_metadata %>% select(Design_Name, Instrument_Model, Library_Name, Library_Layout,
                                                      Library_Source, Library_Selection, Library_Strategy,
                                                      Library_Construction_Protocol) %>% distinct() %>%
                                               arrange(Design_Name)
  # sheet Link_Files_Samples
  dfs_ega[["Link_Files_Samples"]] <- EGA_metadata %>%
    mutate(`Sample alias confidential`=Sample_Id) %>%
    mutate(`Sample alias`=bioSampleId) %>%
    select(Design_Name, `Sample alias`, `Sample alias confidential`, `First Fastq File`, `First Checksum`,
           `First Unencrypted checksum`, `Second Fastq File`, `Second Checksum`, `Second Unencrypted checksum`) %>%
    arrange(Design_Name)

  # sheet DAC
  if (mode=="breast"){
    dfs_ega[["DAC"]] <- bind_rows(tibble(Title="MD, PhD", Name="Fabrice André", Email="fabrice.andre@gustaveroussy.fr",
                                         Telephone="+33 1 42 11 61 59", Organization="Gustave Roussy"),
                                  tibble(Title="MD", Name="Barbara Pistilli", Email="barbara.pistilli@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 23 85", Organization="Gustave Roussy"),
                                  tibble(Title="MD", Name="Fernanda Mosele", Email="fernanda.mosele@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 23 85", Organization="Gustave Roussy"),
                                  tibble(Title="PhD", Name="Marc Deloger", Email="marc.deloger@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 41 26", Organization="Gustave Roussy"),
                                  tibble(Title="PhD", Name="Yoann Pradat", Email="yoann.pradat@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 41 26", Organization="Gustave Roussy"))
  } else if (mode=="lung"){
    dfs_ega[["DAC"]] <- bind_rows(tibble(Title="MD, PhD", Name="Fabrice André", Email="fabrice.andre@gustaveroussy.fr",
                                         Telephone="+33 1 42 11 61 59", Organization="Gustave Roussy"),
                          tibble(Title="MD, PhD", Name="David Planchard", Email="david.planchard@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 56 16", Organization="Gustave Roussy"),
                          tibble(Title="PhD", Name="Guillaume Montagnac", Email="guillaume.montagnac@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 56 16", Organization="Gustave Roussy"),
                                  tibble(Title="PhD", Name="Marc Deloger", Email="marc.deloger@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 41 26", Organization="Gustave Roussy"),
                                  tibble(Title="PhD", Name="Yoann Pradat", Email="yoann.pradat@gustaveroussy.fr",
                                     Telephone="+33 1 42 11 41 26", Organization="Gustave Roussy"))
  }


  dfs_ega
}


main <- function(args){
  outputs_ega_samples <- list(breast=args$output_ega_samples_b,
                              lung=args$output_ega_samples_l)

  outputs_ega_metadata <- list(breast=args$output_ega_metadata_b,
                              lung=args$output_ega_metadata_l)

  for (mode in c("breast", "lung")){
    bio <- load_samples_data(args$cln, args$sam_wes, args$sam_rna, mode)
    stats <- print_statistics(bio)
    ids_anonym <- read.delim(args$ids_anonym, sep="\t")
    bio <- add_anonymous_identifiers(bio, ids_anonym)
    ega_metadata <- make_ega_metadata(bio, stats, mode)

    # # save ega samples
    # cols_sams <- c("bioSampleId", "FASTQ_1", "FASTQ_2")
    # ega_samples <- bio[cols_sams] %>% rename(Sample_Id=bioSampleId)
    # write.table(ega_samples, file=outputs_ega_samples[[mode]], sep="\t", quote=F, row.names=F, na="")
    # cat(paste("-file saved at", outputs_ega_samples[[mode]], "\n"))

    # save ega metadata
    write.xlsx(ega_metadata, file=outputs_ega_metadata[[mode]])
    cat(paste("-file saved at", outputs_ega_metadata[[mode]], "\n"))

    # if some biopsy sites are missing, create sheet to be filled
    if (sum(is.na(bio$Biopsy_Site)) > 0 | sum(bio$Biopsy_Site=="") >0 ){
      fp_mis <- paste0("../../results/tables_paper/ega/ega_samples_", mode, "_biopsy_sites_to_be_filled.xlsx")
      bio_mis <- bio %>% select(Subject_Id, Sample_Id, Sample_Type, Timepoint, Biopsy_Site)
      write.xlsx(bio_mis, fp_mis)
    }
  }

}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Prepare the metadata table for EGA upload.')
  parser$add_argument('--cln', type="character", help='Path to table of clinical data.',
                      default="../../data/cln/curated/cln_icarus_curated_with_qc.tsv")
  parser$add_argument('--sam_wes', type="character", help='Path to table of WES samples and FASTQ filepaths.',
                      default="../../data/wes/config_b1_b2_b3_b4/samples.all.tsv")
  parser$add_argument('--sam_rna', type="character", help='Path to table of RNA-seq samples and FASTQ filepaths.',
                      default="../../data/rna/config/samples.tsv")
  parser$add_argument('--ids_anonym', type="character", help='Path to table of anonymized ids.',
                      default="../../results/tables_paper/table_anonymisation.tsv")
  parser$add_argument('--output_ega_samples_b', type="character", help="Path to output table of EGA samples for breast.",
                      default="../../results/tables_paper/ega/ega_samples_breast.tsv")
  parser$add_argument('--output_ega_samples_l', type="character", help="Path to output table of EGA samples for lung.",
                      default="../../results/tables_paper/ega/ega_samples_lung.tsv")
  parser$add_argument('--output_ega_metadata_b', type="character", help="Path to output table of EGA metadata for breast.",
                      default="../../results/tables_paper/ega/ega_metadata_breast.xlsx")
  parser$add_argument('--output_ega_metadata_l', type="character", help="Path to output table of EGA metadata for lung.",
                      default="../../results/tables_paper/ega/ega_metadata_lung.xlsx")
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)

  main(args)
}

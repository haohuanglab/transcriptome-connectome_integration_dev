library(tidyverse)
library(EWCE) 
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(stringr)
set.seed(1234)

pkg <- tolower("EWCE")
docker_user <- "neurogenomicslab"

flag_rev <- "separate" 

background <- r"[Data/custom_bg.txt]"

target_agepcd <- "3186"
target_dataset <- "8Years"

data_temp <- read_excel("Data/bottom_325genes_tidy.xlsx")
unique_ages <- unique(data_temp$Age)
target_age_label <- unique_ages[grepl(paste0(target_agepcd, "pcd"), unique_ages)]
matching_ages <- data.frame(ages = target_age_label, agepcd = target_agepcd, dataset = target_dataset, genelists = NA)

age <- target_dataset    

columns <- readRDS("Data/celltype_columns.rds")

flag_enrichment <- TRUE

if (flag_enrichment) {
  PLSwhich <- "PLS1" 
  num_gene <- 325     
  reps <- 20000
  metrics = "2Ne_6Lp_7Neloc"
  foldername_celltype <- "lister22.8Years"
  bg_str <- "_background"
  bg_str <- ""

  get_Enrichment <- function(data, bootsrap_output) {
    background <- readLines(background)
      temp <- data[data$Age == matching_ages$ages[1], ]
      g <- temp$`Gene Symbol`
      g <- g[1:num_gene]
      matching_ages$genelists[1] <- paste(g, collapse = ",")
    if (file.exists(bootsrap_output)) {
      bootstrap_list <- readRDS(bootsrap_output)
    } else {
      full_name <- matching_ages$ages[1]
      rda <- paste0("Data/ctd_lister22.", matching_ages$dataset[1], ".rda")
      ctd <- EWCE::load_rdata(rda)
      genelists <- unlist(strsplit(matching_ages$genelists[1], ","))
      full_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd, hits = genelists, 
                                                      bg = background,
                                                      sctSpecies = "human",
                                                      genelistSpecies = "human",
                                                      output_species = "human",
                                                      reps = reps,
                                                      geneSizeControl = FALSE,
                                                      controlledCT = NULL,
                                                      mtc_method = "none",
                                                      annotLevel = 1, verbose = FALSE) 
      bootstrap_list <- list()
      bootstrap_list[[full_name]] <- full_results
      saveRDS(bootstrap_list, file = bootsrap_output)
    }

    result_df <- data.frame(celltype = columns, row.names = columns)
    
    time_point <- names(bootstrap_list)[1]
    cell_types <- bootstrap_list[[time_point]]$results$CellType
    p_values <- bootstrap_list[[time_point]]$results$p
    result_df[as.character(cell_types), as.character(time_point)] <- p_values

    result_df[is.na(result_df)] <- 1
    result_df$celltype <- rownames(result_df)
    return (result_df)
  }

  GRETNA_metric <- substring(metrics, 2, nchar(metrics))  
  folder_celltype <- paste0("Data/", foldername_celltype)
  if (!dir.exists(folder_celltype)) {
    dir.create(folder_celltype, showWarnings = FALSE) 
  }
  
  if (flag_rev == "top") {
    data <- read_excel("Data/top_325genes_tidy.xlsx")
    bootsrap_output <- paste0("Data/full_results_mapping_", "top", num_gene, "genes", bg_str, ".rds")
    enrichment <- get_Enrichment(data, bootsrap_output)
  } else if (flag_rev == "bottom") {
    data <- read_excel("Data/bottom_325genes_tidy.xlsx")
    bootsrap_output <- paste0("Data/full_results_mapping_", "bottom", num_gene, "genes", bg_str, ".rds")
    enrichment <- get_Enrichment(data, bootsrap_output)
  } else if (flag_rev == "separate") {
    data <- read_excel("Data/top_325genes_tidy.xlsx")
    bootsrap_output <- paste0("Data/full_results_mapping_", "top", num_gene, "genes", bg_str, ".rds")
    enrichment_top <- get_Enrichment(data, bootsrap_output)

    data <- read_excel("Data/bottom_325genes_tidy.xlsx")
    bootsrap_output <- paste0("Data/full_results_mapping_", "bottom", num_gene, "genes", bg_str, ".rds")
    enrichment_bottom <- get_Enrichment(data, bootsrap_output)

    enrichment <- rbind(enrichment_top, enrichment_bottom)
  } else {
    message("invalid flag_rev!")
  }

  enrichment_top$gene_set <- "top"
  enrichment_bottom$gene_set <- "bottom"
  
  enrichment <- rbind(enrichment_top, enrichment_bottom)
  enrichment <- enrichment %>% mutate(across(where(is.numeric), p.adjust, method = "fdr"))
  View(enrichment)
  
  if (flag_rev == "both") {
    enrichment_top <- enrichment[1:(nrow(enrichment)/2),]
    enrichment_bottom <- enrichment[(nrow(enrichment)/2 + 1):nrow(enrichment),]
    enrichment_top <- enrichment_top[order(enrichment_top$celltype),]
    enrichment_bottom <- enrichment_bottom[order(enrichment_bottom$celltype),]
    enrichment <- enrichment_top %>% map2_dfr(enrichment_bottom, pmin)
  }
}

write.csv(enrichment, file = paste0("Data/enrichment_top325_", GRETNA_metric, "_", flag_rev, PLSwhich, ".csv"), row.names = FALSE)
write.csv(enrichment_top, file = paste0("Data/enrichment_top.csv"), row.names = FALSE)
write.csv(enrichment_bottom, file = paste0("Data/enrichment_bottom.csv"), row.names = FALSE)

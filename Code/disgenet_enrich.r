library(dplyr)
library(jsonlite)
library(openxlsx)

get.fisher <- function(gene_df1, gene_df2, all_genes, name) {
    result_df <- data.frame(
        subtype1 = character(),
        subtype2 = character(),
        p_value = numeric(),
        odds_ratio = numeric(),
        common_genes = character()
    )
  
    for (subtype1 in names(gene_df1)) {
        for (subtype2 in names(gene_df2)) {
        
        genes1 <- unlist(gene_df1[[subtype1]])
        genes2 <- unlist(gene_df2[[subtype2]])

        genes1 <- intersect(genes1, all_genes)
        genes2 <- intersect(genes2, all_genes)
        
        total_genes <- length(all_genes)
        
        contingency_table <- matrix(0, nrow = 2, ncol = 2)
        
        contingency_table[1, 1] <- length(intersect(genes1, genes2))
        contingency_table[1, 2] <- length(setdiff(genes1, genes2))
        contingency_table[2, 1] <- length(setdiff(genes2, genes1))
        contingency_table[2, 2] <- total_genes - sum(contingency_table)
        
        fisher_result <- fisher.test(contingency_table, alternative = "greater")
      
        result_df <- rbind(result_df, data.frame(
            subtype1 = subtype1,
            subtype2 = subtype2,
            p_value = fisher_result$p.value,
            odds_ratio = fisher_result$estimate,
            common_genes = paste(intersect(genes1, genes2), collapse = ", "),
            percent_common = length(intersect(genes1, genes2))/length(genes1)
        ))
      
    }

  }
  
  return(result_df)
}


background_genes_path <- r"{Data/Homo_sapiens.gene_info}"
background_genes <- read.csv(background_genes_path, sep='\t')
background_genes <- background_genes %>% filter(type_of_gene == "protein-coding")
all_genes <- background_genes$Symbol

disease_gene_path <- r"{Data/disease_genes.json}"
disease_gene_data <- fromJSON(disease_gene_path)

celltype_gene_path <- r"{Data/cell_type_genes.json}"
celltype_gene_data <- fromJSON(celltype_gene_path)

go_gene_path <- r"{Data/genes.json}"
go_gene_data <- fromJSON(go_gene_path)

pls_path_top <- r"{Data/top_325genes_tidy.xlsx}"
pls_path_bottom <- r"{Data/bottom_325genes_tidy.xlsx}"
pls_data_top <- readxl::read_excel(pls_path_top)
pls_data_bottom <- readxl::read_excel(pls_path_bottom)
pls_data <- bind_rows(pls_data_top, pls_data_bottom)


ages <- c(266, 937, 3004, 815, 3369, 6289, 754, 245, 754, 3004)
d0 <- c("ASD", "ASD", "ASD", "BD", "BD", "BD", "MDD", "SCZ", "SCZ", "SCZ")
d1 <- c("Peak 1", "Peak 2", "Peak 3", "Peak 1", "Peak 2", "Peak 3", "Peak 1", "Peak 1", "Peak 2", "Peak 3")
d <- d0

names <- list()
gene_df1 <- list()
for (agepcd in ages) {
    gene_list <- pls_data %>% filter(grepl(paste0(as.character(agepcd),"pcd"), Age)) %>% select(`Gene Symbol`)
    gene_df1[[as.character(agepcd)]] <- gene_list$`Gene Symbol`
    name <- pls_data %>% filter(grepl(paste0(as.character(agepcd),"pcd"), Age)) %>% select(Age)
    name <- gsub(" - .*", "", name)
    name <- gsub("c\\(\"", "", name)
    names[[as.character(agepcd)]] <- name
}

pls_lists <- list()
for (i in 1:length(d)) {
    name <- names[[as.character(ages[i])]]
    name <- gsub(".*\\(", "", name)
    d[i] <- paste0(d[i], " - ", name)
    pls_lists[[d[i]]] <- gene_df1[[as.character(ages[i])]]
}

fisher_dfs <- list()
go_dfs <- list()
for(i in 1:length(d)) {
    disease_subset <- disease_gene_data[grep(d0[i], names(disease_gene_data))]
    names(disease_subset) <- paste0(names(disease_subset), " x ", d1[i])
    disease_subset<- lapply(disease_subset, function(x) {
        x <- intersect(unlist(x), unlist(pls_lists[[d[i]]]))
        return(x)
    })
    fisher_dfs[[d[i]]] <- get.fisher(disease_subset, celltype_gene_data, all_genes, "celltype")
    go_dfs[[d[i]]] <- get.fisher(disease_subset, go_gene_data, all_genes, "go")
}
result <- bind_rows(fisher_dfs)
result2 <- bind_rows(go_dfs)


result$dataset <- gsub(" x .*", "", result$subtype1)
result$dataset <- gsub("_.*", "", result$dataset)
result$subtype1 <- gsub(".*_", "", result$subtype1)

result2$dataset <- gsub(" x .*", "", result2$subtype1)
result2$dataset <- gsub("_.*", "", result2$dataset)
result2$subtype1 <- gsub(".*_", "", result2$subtype1)


result <- result %>% mutate(p_corr = p.adjust(p_value, method = "bonferroni"))
result2 <- result2 %>% mutate(p_corr = p.adjust(p_value, method = "bonferroni"))

result$subtype2 <- gsub("L5-6 NP", "L5-6_NP", result$subtype2)

SIGNIFICANCE <- 1e-3

result_save <- result[c("subtype1", "subtype2", "p_corr", "common_genes")] %>% arrange(subtype1, p_corr)
filtered_data <- result %>% filter(p_corr < SIGNIFICANCE) 


write.xlsx(result_save, "Data/disgenet_results_celltype.xlsx")
write.xlsx(result2,"Data/disgenet_results_bioprocess.xlsx")



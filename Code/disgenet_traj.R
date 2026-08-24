library(dplyr)
library(purrr)
library(ggplot2)
library(jsonlite)
library(ggpubr)
library(xlsx)
library(reshape2)
library(tidyr)
library(latex2exp)
library(stringr)

standardize.list <- function(genelist, all_genes) {
  
  synonym_list <- all_genes$synonym
  bg_list <- all_genes$Symbol
  all_genes_list <- union(synonym_list, bg_list)
  genelist <- intersect(genelist, all_genes_list)
  for (i in 1:length(genelist)) {
    if (genelist[i] %in% synonym_list) {
      genelist[i] <- bg_list[which(synonym_list == genelist[i])]
    }
  }
  return (genelist)
  
}

get.fisher <- function(genelist1, geenelist2, age, disease, all_genes) {
  
  result_df <- data.frame(
    age = character(),
    disease = character(), 
    p_value = numeric(),
    odds_ratio = numeric(),
    common_genes = character()
  )
  
  genes1 <- intersect(genelist1, all_genes)
  genes2 <- intersect(genelist2, all_genes)
  
  total_genes <- length(all_genes)
  
  contingency_table <- matrix(0, nrow = 2, ncol = 2)
  
  contingency_table[1, 1] <- length(intersect(genes1, genes2))
  contingency_table[1, 2] <- length(setdiff(genes1, genes2))
  contingency_table[2, 1] <- length(setdiff(genes2, genes1))
  contingency_table[2, 2] <- total_genes - sum(contingency_table)
  
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  return(data.frame(
    age = age,
    disease = disease,
    p_value = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    common_genes = paste(intersect(genes1, genes2), collapse = ", "),
    percent_common = length(intersect(genes1, genes2))/length(genes1)))
}

disease_genes_path <- r"{Data\Disgenet_genes.xlsx}"
sheet_names <- c("C1510586", "C0005586", "C1269683", "C0036341")
names <- c("ASD", "BD", "MDD", "SCZ")

disease_genes <- lapply(sheet_names, function(sheet_name) { 
  read.xlsx(disease_genes_path, sheet_name)$Gene
})
names(disease_genes) <- names

pls_path <- "\Data\both_650genes_tidy.xlsx"
pls_data <- readxl::read_excel(pls_path)
gene_df1 <- list()

for (age in unique(pls_data$Age)) {
  gene_list <- pls_data %>% filter(Age == age) %>% select(`Gene Symbol`) %>% pull()
  agepcd <- gsub(".* - ", "", age)
  gene_df1[[agepcd]] <- gene_list
}


background_genes_path <- r"[\Data\Homo_sapiens.gene_info]"
background_genes <- read.csv(background_genes_path, sep='\t')

all_genes <- background_genes %>% filter(type_of_gene == "protein-coding")
all_genes <- all_genes$Symbol


disease_df <- data.frame()
for (disease in names(disease_genes)) {
  genelist2 <- disease_genes[[disease]]
  for (age in names(gene_df1)) {
    genelist1 <- gene_df1[[age]]
    result_df <- get.fisher(genelist1, genelist2, age, disease, all_genes)
    disease_df <- rbind(disease_df, result_df)
  }
}

disease_df$p_correct <- p.adjust(disease_df$p_value, method = "bonferroni")
disease_df$agepcd <- gsub(".* - ", "", disease_df$age)
disease_df$agepcd <- gsub("pcd.*", "", disease_df$agepcd)
disease_df$agepcd <- as.numeric(disease_df$agepcd)
disease_df <- disease_df %>% mutate(p_correct = ifelse(p_correct < 1e-12, 1e-12, p_correct))
disease_df$significant <- disease_df$p_correct < 1e-5

folder_main_geneWeights <- "Data"  

filename_genesraw_list = c("ASD", "BD", "MDD", "SCZ")
num_gene <- 650

df_xtick_full <- data.frame(agepcd = 38*7 + c(0,1,2,3,4,6,8,15,22)*365,
                            tick = c("Birth","1","2","3", "4", "6", "8","15 ", "22"))
max_y_value <- 12

k_value = 30
bs_value = "cs"
span = 0.2

disease_df$pvalue <- disease_df$p_value
disease_df$padj <- disease_df$p_correct
disease_df$method <- "disgenet"


disease_df_copy <- disease_df[,c("agepcd","padj","disease", "pvalue", "method")]

max_y_value <- 12

h_ggplot <- ggplot(data = disease_df, aes(x = agepcd, y = -log10(p_correct))) +
  geom_point(aes(color=significant), size=3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  geom_smooth(data=disease_df, aes(x = agepcd, y = -log10(p_correct)), method="gam", se=FALSE, color="#FF6961", formula = y ~ s(log10(x), k = k_value, bs=bs_value), n=1000, alpha=0.7) +
  coord_trans(x = "log10", xlim = c(min(df_xtick_full$agepcd)-14, max(df_xtick_full$agepcd)+28)) +
  scale_x_continuous(breaks = df_xtick_full$agepcd, labels = df_xtick_full$tick) + 
  scale_y_continuous(limits = c(0, max_y_value), breaks = seq(0, max_y_value, length.out = 7)) +  
  facet_wrap(disease ~ ., scales = "free_y", ncol = 1) +
  labs(title = NULL,
       x = "Age",
       y = TeX("$-log_{10} pBonferroni$")) +
  theme_bw() +
  theme(text = element_text(size = 20),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20),  
        axis.text.x = element_text(size = 18),  
        axis.text.y = element_text(size = 18),  
        strip.text = element_text(size = 18),  
        strip.background = element_blank(),  
        legend.position = "none")

max_y_value <- 5

write.xlsx(disease_df, "Data\disgenet_traj_results.xlsx")


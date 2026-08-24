if (!require("readxl")) install.packages("readxl")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("gridExtra")) install.packages("gridExtra")

library(readxl)
library(tidyverse)
library(gridExtra)

genes_of_interest <- c(
  "HTR1A", "HTR1B", "HTR2A", "HTR4", "HTR6", "SLC6A4",
  "DRD1", "DRD2", "SLC6A3",
  "GABRA1", "GABRB2", "GABRG2",
  "CHRNA4", "CHRNB2", "CHRM1",
  "HRH3",
  "GRIN1", "GRIN2A", "GRIN2B",
  "SLC6A2", "SLC18A3",
  "GRM5",
  "CNR1", "OPRM1"
)


gene_colors <- c(
  # Serotonin (reds)
  "HTR1A"  = "#990000",  # deep red
  "HTR1B"  = "#cc3300",  # strong brick
  "HTR2A"  = "#e6550d",  # orange-red
  "HTR4"   = "#fd8d3c",  # tangerine
  "HTR6"   = "#fdae6b",  # warm peach
  "SLC6A4" = "#fdd0a2",  # soft peach

  # Dopamine (blues)
  "DRD1"   = "#08306b",  # very dark blue
  "DRD2"   = "#2171b5",  # medium blue
  "SLC6A3" = "#6baed6",  # light sky blue

  # GABA (greens)
  "GABRA1" = "#00441b",  # forest green
  "GABRB2" = "#238b45",  # medium green
  "GABRG2" = "#74c476",  # minty green

  # Acetylcholine (purples/pinks)
  "CHRNA4" = "#49006a",  # very dark purple
  "CHRNB2" = "#7a0177",  # deep purple
  "CHRM1"  = "#c51b8a",  # strong magenta

  # Histamine (yellows)
  "HRH3"   = "#ffd700",  # gold/yellow

  # NMDA (oranges)
  "GRIN1"  = "#f16913",  # strong orange
  "GRIN2A" = "#fdae6b",  # peach-orange
  "GRIN2B" = "#fdd0a2",  # light peach

  # Glutamate (greys)
  "GRM5"   = "#252525",  # dark grey

  # Noradrenergic (violet shades)
  "SLC6A2"  = "#3f007d",  # indigo
  "SLC18A3" = "#807dba",  # soft violet

  # Cannabinoid & Opioid (browns)
  "CNR1"   = "#8c510a",   # dark ochre
  "OPRM1"  = "#d8b365"    # tan
)


df_xtick <- data.frame(
  age = 38 * 7 + c(0,1,2,3,4,6,8,15,22)*365,
  tick = c("Birth", "1", "2", "3", "4", "6", "8", "15  ", "  22")
)

extract_gene_ranks <- function(file_path, genes, max_rank = 326) {
  raw_data <- read_excel(file_path, col_names = FALSE)
  rank_df <- data.frame()

  for (i in seq(1, ncol(raw_data), by = 3)) {
    age_label <- as.character(raw_data[1, i])
    
    age_string <- str_trim(str_extract(age_label, "(?<=-).*"))
    age_pcd <- as.numeric(str_extract(age_string, "\\d+"))

    gene_symbols <- raw_data[-1, i + 1] %>% unlist() %>% as.character()
    ranks <- setNames(seq_along(gene_symbols), gene_symbols)

    age_ranks <- tibble(
      Gene = genes,
      Rank = ifelse(genes %in% names(ranks), ranks[genes], max_rank),
      AgeLabel = age_label,
      AgeDays = age_pcd
    )

    rank_df <- bind_rows(rank_df, age_ranks)
  }
  return(rank_df)
}

library(ggrepel)  

library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

plot_gene_ranks <- function(rank_df, title_text) {
    plot_data <- rank_df

    plot_data$Gene <- factor(
        str_trim(plot_data$Gene),
        levels = genes_of_interest  
        )
  
  peak_labels <- plot_data %>%
    group_by(Gene) %>%
    filter(Rank == min(Rank)) %>%
    ungroup() %>%
    filter(Rank < 325)
  peak_labels <- peak_labels %>% filter(!is.na(Gene), Gene != "")

    ggplot(plot_data, aes(x = AgeDays, y = Rank, color = Gene, group = Gene)) +
    geom_line(size = 1) +
    geom_point() +
    scale_x_log10(
      breaks = df_xtick$age,
      labels = df_xtick$tick
    ) +
    scale_y_reverse(
      breaks = c(1, 50, 100, 150, 200, 250, 300, 326),
      labels = c("1", "50", "100", "150", "200", "250", "300", ">325")
    ) +
    scale_color_manual(
      values = gene_colors,
      breaks = levels(plot_data$Gene)  
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(margin = margin(b = 4)),
      legend.key.height = unit(1.0, "lines"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.1, "cm"),
      panel.grid.major.y = element_line(color = "grey80"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey80"),
      panel.grid.minor.x = element_blank()
    ) +
    guides(color = guide_legend(ncol = 1)) +
    labs(
      title = title_text,
      x = "Age (years)",
      y = "PLS Ranking (1 = Higher Rank)",
      color = "Gene"
    ) + geom_label_repel(
        data = peak_labels,
        aes(label = as.character(Gene)), 
        size = 3,
        fill = alpha("white", 0.7),
        label.size = NA,
        label.r = unit(0, "pt"),
        box.padding = 0.3,
        segment.color = NA,
        min.segment.length = 0,
        show.legend = FALSE
    ) 
}


top_file <- "Data\top_325genes.xlsx"
bottom_file <- "Data\bottom_325genes.xlsx"

rank_df_top <- extract_gene_ranks(top_file, genes_of_interest)
rank_df_bottom <- extract_gene_ranks(bottom_file, genes_of_interest)

p1 <- plot_gene_ranks(rank_df_top, "Top 325 PLS Neurotransmitter Gene Rankings Across Ages")
p2 <- plot_gene_ranks(rank_df_bottom, "Bottom 325 PLS Neurotransmitter Gene Rankings Across Ages")

if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
library(cowplot)

p2_no_legend <- p2 + theme(legend.position = "none")

legend <- get_legend(p1)

combined_plot <- plot_grid(
  plot_grid(p1 + theme(legend.position = "none"), p2_no_legend, ncol = 1, align = "v"),
  legend,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

print(combined_plot)
ggsave(
  filename = "Data\neurotransmitter_gene_rankings.pdf",
  plot = combined_plot,
  width = 8,     
  height = 6, 
  dpi = 900   
)

find_always_low_rank_genes <- function(rank_df, max_rank = 325) {
  always_low <- rank_df %>%
    group_by(Gene) %>%
    summarize(all_low = all(Rank > max_rank)) %>%
    filter(all_low) %>%
    pull(Gene)
  
  return(always_low)
}

always_low_top <- find_always_low_rank_genes(rank_df_top)
always_low_bottom <- find_always_low_rank_genes(rank_df_bottom)

always_low_both <- intersect(always_low_top, always_low_bottom)

cat("Genes always ranked > 325 in BOTH datasets:\n")
print(always_low_both)



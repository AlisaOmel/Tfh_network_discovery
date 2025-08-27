library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)
library(msigdbr)
library(effsize)  # for Cliff's delta
library(pheatmap)
library(RColorBrewer)
library(tibble)
library(ggrepel)

#Load in the cells
###########################################################################
combined <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_nonTfh_Tfh_combo.rds")
hiv.groups <- c('PLWHIV')
combined <- subset(combined, subset = HIV %in% hiv.groups)
###########################################################################

#Load msigDB pathway 
###########################################################################
GS.hallmark <- getGeneSets(library = "H")

###########################################################################

#Run analysis
###########################################################################
ssg_seurat = runEscape(combined, method = "ssGSEA", gene.sets = GS.hallmark, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)

# Reorder the levels of the cell_type metadata
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("non-Tfh", "Tfh"))

ssg_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/scGSEA_hallmark_pathway_analysis_CM_PLWHIV.RSD")
###########################################################################



#Loop through all pathways and make a summary statistic and saves
########################################################################################################
# Extract enrichment matrix
enrichment_matrix <- as.data.frame(t(ssg_seurat[["escape.ssGSEA"]]@data))

# Add metadata for cell type classification
enrichment_matrix$cell_type <- ssg_seurat$cell_type

# Initialize results list
results <- list()

# Loop through each pathway
for (pathway in colnames(enrichment_matrix)[!colnames(enrichment_matrix) %in% "cell_type"]) {
  
  # Subset scores by cell type
  tfh_scores <- enrichment_matrix %>% filter(cell_type == "Tfh") %>% pull(pathway)
  non_tfh_scores <- enrichment_matrix %>% filter(cell_type == "non-Tfh") %>% pull(pathway)
  
  # Compute means
  mean_tfh <- mean(tfh_scores, na.rm = TRUE)
  mean_non_tfh <- mean(non_tfh_scores, na.rm = TRUE)
  
  # Cliff's delta
  cliff <- tryCatch({
    effsize::cliff.delta(tfh_scores, non_tfh_scores)$estimate
  }, error = function(e) NA)
  
  # Wilcoxon p-value
  pval <- tryCatch({
    wilcox.test(tfh_scores, non_tfh_scores)$p.value
  }, error = function(e) NA)
  
  # Append to results
  results[[pathway]] <- data.frame(
    Pathway = pathway,
    Mean_Tfh = mean_tfh,
    mean_non_tfh = mean_non_tfh,
    Cliffs_Delta = cliff,
    P_Value = pval
  )
}

# Combine and write CSV
summary_df <- do.call(rbind, results)

print(summary_df)

write.csv(summary_df, "/ix3/djishnu/Alisa/Tfh/scRNA/hallmark_ssGSEA_GCTfh_vs_nonTfh_summary_CM_PLWHIV_TEST.csv", row.names = FALSE)


#Heat Plot
########################################################################
ssg_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/scGSEA_PID_pathway_analysis_Th17_PLWHIV.RSD")
summary_df <- read.csv("/ix3/djishnu/Alisa/Tfh/scRNA/pid_ssGSEA_GCTfh_vs_nonTfh_summary_Th17_PLWHIV.csv")

# Get top 15 by absolute value
top15 <- summary_df %>%
  arrange(desc(abs(Cliffs_Delta))) %>%
  slice(1:15) %>%
  mutate(Pathway = factor(Pathway, levels = rev(Pathway)))  # order from largest to smallest

# Add a dummy column for the x-axis (single column)
top15$column <- "Cliff's Delta"

ggplot(top15, aes(x = column, y = Pathway, fill = Cliffs_Delta)) +
  geom_tile(color = "white", width = 0.3, height = 0.8) +  # narrower width
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0, name = "Cliff's Δ"
  ) +
  theme_minimal() +
  labs(
    #title = "Top 15 Pathways by Absolute Cliff’s Delta",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )

########################################################################
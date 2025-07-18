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

#Load scGSEA results
###########################################################################################
#ssg_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/scGSEA_PID_pathway_analysis_Th17_PLWHIV.RSD")
#summary_df <- read.csv("/ix3/djishnu/Alisa/Tfh/scRNA/pid_ssGSEA_GCTfh_vs_nonTfh_summary_Th17_PLWHIV.csv")

#ssg_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/scGSEA_hallmark_pathway_analysis_Th17_PLWHIV.RSD")
#summary_df <- read.csv("/ix3/djishnu/Alisa/Tfh/scRNA/hallmark_ssGSEA_GCTfh_vs_nonTfh_summary_Th17_PLWHIV.csv")

ssg_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/scGSEA_KEGG_pathway_analysis_Th17_PLWHIV.RSD")
summary_df <- read.csv("/ix3/djishnu/Alisa/Tfh/scRNA/KEGG_ssGSEA_GCTfh_vs_nonTfh_summary_Th17_PLWIH.csv")
###########################################################################################

#Plot
###########################################################################################
# Specify gene set to plot
#gene_name <- "PID_IL12_2PATHWAY"
#gene_name <- "PID_IL27_PATHWAY"
#gene_name <- "PID_IL23_PATHWAY"
#gene_name <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
gene_name <- "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"

converted_name <- gsub("_", "-", gene_name)

# Get ssGSEA enrichment matrix
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

# Check if gene set exists
if (!(converted_name %in% rownames(ssgsea_mat))) {
  stop("Gene set not found in ssGSEA matrix.")
}

# Create enrichment score data frame
df <- data.frame(
  score = ssgsea_mat[converted_name, ],
  group = ssg_seurat$cell_type
)
df <- df[df$group %in% c("non-Tfh", "Tfh"), ]
df$group <- factor(df$group, levels = c("non-Tfh", "Tfh"))

# Stop if no variability
if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
  stop("No variability in enrichment scores.")
}

# Calculate Cliff's Delta
cliff <- tryCatch({
  cidv2(
    df[df$group == "non-Tfh", "score"],
    df[df$group == "Tfh", "score"],
    alpha = 0.05
  )[c("d.hat", "p.value")]
}, error = function(e) {
  stop("Cliff's Delta calculation failed.")
})

# Calculate group means
avg_non_tfh <- mean(df$score[df$group == "non-Tfh"], na.rm = TRUE)
avg_tfh <- mean(df$score[df$group == "Tfh"], na.rm = TRUE)

# Prepare plotting data
df_plot <- as.data.frame(t(as.matrix(ssg_seurat@assays$escape.ssGSEA@data)))
df_plot$cell_type <- ssg_seurat$cell_type
df_plot <- df_plot[, c(converted_name, "cell_type")]
colnames(df_plot)[1] <- "score"

df_for_points <- df_plot %>%
  group_by(cell_type) %>%
  sample_frac(0.7) %>%
  ungroup()

df_plot$cell_type <- factor(df_plot$cell_type, levels = c("Tfh", "non-Tfh"))
levels(df_plot$cell_type)[levels(df_plot$cell_type) == "non-Tfh"] <- "Th17"
df_for_points$cell_type <- factor(df_for_points$cell_type, levels = c("Tfh", "non-Tfh"))
levels(df_for_points$cell_type)[levels(df_for_points$cell_type) == "non-Tfh"] <- "Th17"

# Plot label
delta_label <- bquote("Cliff's " * Delta * " = " * .(round(cliff[["d.hat"]], 2)))

my_comparison <- list(c("Th17", "Tfh"))

# Create plot
gs <- ggplot(df_plot, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
  geom_jitter(data = df_for_points, width = 0.2, size = 1, alpha = 0.1) +
  theme_classic() +
  labs(x = "Cell Type", y = paste0(gene_name, " Enrichment Score")) +
  stat_compare_means(comparisons = my_comparison, 
                     method = "wilcox.test", 
                     label = "p.signif",
                     size = 6.5) +
  annotate("text", x = 1.5, y = max(df_plot$score, na.rm = TRUE), 
           label = delta_label, vjust = 0.5, size = 6.5) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  )

# Save to file
ggsave(
  filename = paste0("/ix3/djishnu/Alisa/Tfh/scRNA/modules_gs/final_plots/gut_geyser_Th17_PLWHIV_", gene_name, "_boxplot_smaller.pdf"),
  plot = gs,
  device = "pdf",
  width = 7,
  height = 5
)

# Optional: Print the plot
print(gs)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)
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

#Use this to match your gene lists to make sure the gene is in your set before running enrichment
seurat_genes <- rownames(combined)
###########################################################################

#Load in the sets of interest - I had a bunch of different .csv, .gmt, .txt etc
###########################################################################

#Modules
pps_df  <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/process_network_prop_output/PPS_1_significant_gene_list.csv')
pps <- as.list(pps_df$'Genes')
pps_genes <- intersect(unlist(pps), seurat_genes)


taiji_df  <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/processed_Taiji/final_genes_all_sets_152.csv')
taiji <- as.list(taiji_df$'Genes')
taiji_genes <- intersect(unlist(taiji), seurat_genes)

vinuesa_df  <-read.csv("/ix/djishnu/Alisa/Tfh/vinuesa_list.csv")
vinuesa <- as.list(vinuesa_df$'Genes')
vinuesa_genes <- intersect(unlist(vinuesa), seurat_genes)

tfh_review_df  <- read.csv("/ix/djishnu/Alisa/Tfh/TFH_review_list.csv")
tfh_review <- as.list(tfh_review_df$'Genes')
tfh_review_genes <- intersect(unlist(tfh_review), seurat_genes)
###########################################################################


#Make a list of all gene lists of interest
###########################################################################
pathway_names <- c("PPI", "GRN", "Review-2022", "Review-2016")
gene_set <- list(
  PPI = pps_genes,
  GRN = taiji_genes,
  Review_2022 = tfh_review_genes,
  Review_2016 = vinuesa_genes
)
###########################################################################


#Run ssGSEA
###########################################################################
ssg_seurat = runEscape(combined, method = "ssGSEA", gene.sets = gene_set, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)


DefaultAssay(ssg_seurat)  <- "escape.ssGSEA"

# Reorder in the way you want it to show up on your plot
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("non-Tfh", "Tfh"))
saveRDS(ssg_seurat, file = "/ix3/djishnu/Alisa/Tfh/scGSEA_modules_pathway_analysis_CM_PLWHIV.RSD")

#ssg_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/scGSEA_modules_pathway_analysis_CM_PLWHIV.RSD")

#Gives the comparison for the p-value calculation in the plot
my_comparison <- list(c("CM non-Tfh CD4", "Tfh"))
# Convert underscores to hyphens to match ssGSEA output (Anything with a hyphen in the pathway gets changed by ssGSEA from _ to -)
valid_gene_sets <- names(gene_set)
converted_names <- gsub("_", "-", valid_gene_sets)

# Check which ones exist in the ssGSEA data
existing_sets <- rownames(ssg_seurat[["escape.ssGSEA"]]@data)
valid_gene_sets <- valid_gene_sets[converted_names %in% existing_sets]
###########################################################################


#Loop through gene sets and plot each one saving a summary.csv
###########################################################################
# Prepare results list
results_list <- list()

# Get ssGSEA enrichment matrix once
ssgsea_mat <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")

for (gene_name in valid_gene_sets) {
  cat("Plotting:", gene_name, "\n")
  
  converted_name <- gsub("_", "-", gene_name)
  
  # Check if gene set is present in assay
  if (!(converted_name %in% rownames(ssgsea_mat))) {
    cat("  Skipping - gene set not found\n")
    next
  }
  
  # Build data frame of enrichment scores and group labels
  df <- data.frame(
    score = ssgsea_mat[converted_name, ],
    group = ssg_seurat$cell_type
  )
  df <- df[df$group %in% c("non-Tfh", "Tfh"), ]
  df$group <- factor(df$group, levels = c("non-Tfh", "Tfh"))
  
  # Skip if invalid data
  if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
    cat("  Skipping - no variability\n")
    next
  }
  
  # Run Cliff's Delta
  cliff <- tryCatch({
    cidv2(
      df[df$group == "non-Tfh", "score"],
      df[df$group == "Tfh", "score"],
      alpha = 0.05
    )[c("d.hat", "p.value")]
  }, error = function(e) {
    cat("  Skipping - cidv2 failed\n")
    return(NULL)
  })
  
  if (is.null(cliff)) next
  
  # Compute group means
  avg_th17 <- mean(df$score[df$group == "non-Tfh"], na.rm = TRUE)
  avg_tfh <- mean(df$score[df$group == "Tfh"], na.rm = TRUE)
  
  # Save result for summary
  results_list[[gene_name]] <- data.frame(
    gene_set = gene_name,
    d.hat = cliff[["d.hat"]],
    p.value = cliff[["p.value"]],
    mean_Th17 = avg_th17,
    mean_Tfh = avg_tfh
  )
  
  # Plot label
  delta_label <- bquote("Cliff's " * Delta * " = " * .(round(cliff[["d.hat"]], 2)))
  
  #Make Boxplots: # Extract scores from the assay
  df <- as.data.frame(t(as.matrix(ssg_seurat@assays$escape.ssGSEA@data)))
  df$cell_type <- ssg_seurat$cell_type
  df <- df[, c(converted_name, "cell_type")]
  colnames(df)[1] <- "score"
  
  df_for_points <- df %>%
    group_by(cell_type) %>%
    sample_frac(0.7) %>%  # Sample 10% of each group
    ungroup()
  
  df$cell_type <- factor(df$cell_type, levels = c("Tfh", "non-Tfh"))
  levels(df$cell_type)[levels(df$cell_type) == "non-Tfh"] <- "CM non-Tfh CD4"
  
  df_for_points$cell_type <- factor(df_for_points$cell_type, levels = c("Tfh", "non-Tfh"))
  levels(df_for_points$cell_type)[levels(df_for_points$cell_type) == "non-Tfh"] <- "CM non-Tfh CD4"
  
  gs <- ggplot(df, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8) +
    geom_jitter(data = df_for_points, width = 0.2, size = 1, alpha = 0.1) +
    theme_classic() +
    labs(x = "Cell Type", y = paste0(gene_name, " Enrichment Score")) +
    stat_compare_means(comparisons = my_comparison, 
                       method = "wilcox.test", 
                       label = "p.signif",
                       size=6.5) +
    annotate("text", x = 1.5, y = max(df$score, na.rm = TRUE), 
             label = delta_label, vjust = 0.5, size = 6.5) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  gs
  
  # Save plot
  ggsave(filename = paste0("/ix3/djishnu/Alisa/Tfh/scRNA/modules_gs/gut_geyser_CM_PLWHIV_", gene_name, "_boxplot_smaller.pdf"),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}

# Write summary CSV
if (length(results_list) > 0) {
  summary_df <- do.call(rbind, results_list)
  write.csv(summary_df,
            "/ix3/djishnu/Alisa/Tfh/scRNA/modules_gs/cliffs_delta_summary_Tfh_vs_CM_PLWHIV.csv",
            row.names = FALSE)
  cat("Summary CSV saved.\n")
} else {
  cat("No valid results to write to summary.\n")
}

##################################################################################





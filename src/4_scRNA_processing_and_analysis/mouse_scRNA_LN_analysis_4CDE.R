library(zellkonverter)
library(Seurat)
library(SummarizedExperiment)
library(ggplot2)

library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(BiocParallel)
library(rogme)

#Read and convert mouse scRNA dataset
###########################################################################
sce <- readH5AD("/ix3/djishnu/Alisa/Tfh/mouse_scRNA/zhongli_ref_202401203_mannually_woDoublet.h5ad")

assayNames(sce)   
colnames(colData(sce))

table(colData(sce)$cell_type)
assayNames(sce)


seurat_obj <- as.Seurat(sce, counts = "counts", data = "X")
DefaultAssay(seurat_obj) <- "RNA"  # this will now point to the correct data

seurat_obj <- as.Seurat(sce, counts = "counts", data = "X")
saveRDS(seurat_obj, file = "/ix3/djishnu/Alisa/Tfh/mouse_scRNA/mouse_scRNA.rds")

###########################################################################


#Subset to groups of interest
###########################################################################
seurat_obj <- readRDS('/ix3/djishnu/Alisa/Tfh/mouse_scRNA/mouse_scRNA.rds')

subset_seurat <- subset(seurat_obj, subset = cell_type %in% c("Tfh", "Resting T"))


subset_seurat[["RNA"]] <- subset_seurat[["originalexp"]]
#subset_seurat[["originalexp"]] <- NULL
DefaultAssay(subset_seurat) <- "RNA"


#Use this to match your gene lists to make sure the gene is in your set before running enrichment
seurat_genes <- rownames(subset_seurat)

DimPlot(subset_seurat, group.by = "cell_type", label = TRUE) +
  ggtitle("UMAP - All Cell Types")


Assays(subset_seurat)
DefaultAssay(subset_seurat) <- "originalexp"

subset_seurat[["originalexp"]]@data[1:5, 1:5]

#saveRDS(ssg_seurat, file = "/ix3/djishnu/Alisa/Tfh/mouse_scRNA/scGSEA_modules_pathway_analysis_RestingT.RSD")
#subset_seurat <- readRDS("/ix3/djishnu/Alisa/Tfh/mouse_scRNA/scGSEA_modules_pathway_analysis_RestingT.RSD")

#Load in the sets of interest - I had a bunch of different .csv, .gmt, .txt etc
###########################################################################


format_string <- function(strings) {
  sapply(strings, function(x) {
    # Find position of first digit
    digit_pos <- regexpr("\\d", x)[1]
    
    if (digit_pos == -1) {
      # No number found: keep first letter capital, rest lowercase
      return(paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))))
    }
    
    first_part <- substr(x, 1, digit_pos - 1)
    digit <- substr(x, digit_pos, digit_pos)
    rest <- substr(x, digit_pos + 1, nchar(x))
    
    formatted <- paste0(
      toupper(substr(first_part, 1, 1)),  # keep first letter capital
      tolower(substr(first_part, 2, nchar(first_part))),  # lowercase rest of prefix
      digit,
      tolower(rest)
    )
    return(formatted)
  })
}



crotty_df  <-read.csv("/ix/djishnu/Alisa/Tfh/croty_list.csv")
crotty <- as.list(crotty_df$'Genes')
crotty_format <- format_string(crotty)
crotty_genes <- intersect(unlist(crotty_format), seurat_genes)

vinuesa_df  <-read.csv("/ix/djishnu/Alisa/Tfh/vinuesa_list.csv")
vinuesa <- as.list(vinuesa_df$'Genes')
vinuesa_format <- format_string(vinuesa)
vinuesa_genes <- intersect(unlist(vinuesa_format), seurat_genes)

tfh_review_df  <- read.csv("/ix/djishnu/Alisa/Tfh/TFH_review_list.csv")
tfh_review <- as.list(tfh_review_df$'Genes')
tfh_review_format <- format_string(tfh_review)
tfh_review_genes <- intersect(unlist(tfh_review_format), seurat_genes)

bcell_df  <- read.csv("/ix/djishnu/Alisa/Tfh/beckys_list.csv")
bcell <- as.list(bcell_df$'Genes')
bcell_genes_format <- format_string(bcell)
bcell_genes <- intersect(unlist(bcell_genes_format), seurat_genes)

pps_df  <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/process_network_prop_output/PPS_1_significant_gene_list.csv')
pps <- as.list(pps_df$'Genes')
pps_format <- format_string(pps)
pps_genes <- intersect(unlist(pps_format), seurat_genes)


taiji_df  <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/processed_Taiji/final_genes_all_sets_152.csv')
taiji <- as.list(taiji_df$'Genes')
taiji_format <- format_string(taiji)
taiji_genes <- intersect(unlist(taiji_format), seurat_genes)

###########################################################################

#Make a list of all gene lists of interest
###########################################################################
pathway_names <- c("Review-2014", "Bcell", "PPI", "Taiji", "Review-2016", "Review-2022")
gene_set <- list(
  Review_2014 = crotty_genes,
  Bcell = bcell_genes,
  PPI = pps_genes,
  Taiji = taiji_genes,
  Review_2016 = vinuesa_genes,
  Review_2022 = tfh_review_genes
)
###########################################################################


#Run ssGSEA
###########################################################################
ssg_seurat = runEscape(subset_seurat, method = "ssGSEA", gene.sets = gene_set, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)


DefaultAssay(ssg_seurat)  <- "escape.ssGSEA"

# Reorder in the way you want it to show up on your plot
ssg_seurat$cell_type <- factor(ssg_seurat$cell_type, levels = c("Tfh", "Resting T"))

saveRDS(ssg_seurat, file = "/ix3/djishnu/Alisa/Tfh/mouse_scRNA/scGSEA_modules_pathway_analysis_mouse_Tfh_RestingT.RSD")

#Gives the comparison for the p-value calculation in the plot
my_comparison <- list(c("Tfh", "Resting T"))


# Convert underscores to hyphens to match ssGSEA output (Anything with a hyphen in the pathway gets changed by ssGSEA from _ to -)
valid_gene_sets <- names(gene_set)
converted_names <- gsub("_", "-", valid_gene_sets)

# Check which ones exist in the ssGSEA data
existing_sets <- rownames(ssg_seurat[["escape.ssGSEA"]]@data)
valid_gene_sets <- valid_gene_sets[converted_names %in% existing_sets]
###########################################################################




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
  df <- df[df$group %in% c("Resting T", "Tfh"), ]
  df$group <- factor(df$group, levels = c("Resting T", "Tfh"))
  
  # Skip if invalid data
  if (all(is.na(df$score)) || length(unique(df$score)) <= 1) {
    cat("  Skipping - no variability\n")
    next
  }
  
  # Run Cliff's Delta
  cliff <- tryCatch({
    cidv2(
      df[df$group == "Resting T", "score"],
      df[df$group == "Tfh", "score"],
      alpha = 0.05
    )[c("d.hat", "p.value")]
  }, error = function(e) {
    cat("  Skipping - cidv2 failed\n")
    return(NULL)
  })
  
  if (is.null(cliff)) next
  
  # Compute group means
  avg_th17 <- mean(df$score[df$group == "Resting T"], na.rm = TRUE)
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
  delta_label <- bquote("Cliff's " * Delta * " = " * .(-1*round(cliff[["d.hat"]], 2)))
  
  df <- as.data.frame(t(as.matrix(ssg_seurat@assays$escape.ssGSEA@data)))
  df$cell_type <- ssg_seurat$cell_type
  df <- df[, c(converted_name, "cell_type")]
  colnames(df)[1] <- "score"
  
  df_for_points <- df %>%
    group_by(cell_type) %>%
    sample_frac(1.0) %>%  # Sample 10% of each group
    ungroup()
  
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
  ggsave(filename = paste0("/ix3/djishnu/Alisa/Tfh/mouse_scRNA/final_plots/scRNA/mouse_geyser_Tfh_vs_RestingT", gene_name, "boxplot_smaller.pdf"),
         plot = gs,
         device = "pdf",
         width = 7,
         height = 5)
}

# Write summary CSV
if (length(results_list) > 0) {
  summary_df <- do.call(rbind, results_list)
  write.csv(summary_df,
            "/ix3/djishnu/Alisa/Tfh/mouse_scRNA/cliffs_delta_summary_Tfh_vs_RestingT.csv",
            row.names = FALSE)
  cat("Summary CSV saved.\n")
} else {
  cat("No valid results to write to summary.\n")
}

##################################################################################

library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)
library(escape)
library(ggplot2)
library(BiocParallel)
library(rogme)
library(Cairo)

seurat_obj <- readRDS('/ix3/djishnu/Alisa/Tfh/scRNA/spatial/mouse_scRNA_spatial_koushul.rds')


DimPlot(seurat_obj, reduction = 'spatial', group.by = "cell_type", label = F) #+
 # ggtitle("UMAP - All Cell Types")

subset_seurat <- subset(seurat_obj, subset = cell_type %in% c("T_follicular_helper", "T_CD4"))


seurat_obj[["RNA"]] <- seurat_obj[["originalexp"]]
DefaultAssay(seurat_obj) <- "RNA"

#seurat_genes <- rownames(seurat_obj)
seurat_genes <- rownames(subset_seurat)

DefaultAssay(seurat_obj) <- "originalexp"

seurat_obj[["originalexp"]]@data[1:5, 1:5]

#Load in the sets of interest - I had a bunch of different .csv, .gmt, .txt etc
###########################################################################
vinuesa_df  <-read.csv("/ix/djishnu/Alisa/Tfh/vinuesa_list.csv")
vinuesa <- as.list(vinuesa_df$'Genes')
vinuesa_genes <- intersect(unlist(vinuesa), seurat_genes)

tfh_review_df  <- read.csv("/ix/djishnu/Alisa/Tfh/TFH_review_list.csv")
tfh_review <- as.list(tfh_review_df$'Genes')
tfh_review_genes <- intersect(unlist(tfh_review), seurat_genes)

pps_df  <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/process_network_prop_output/PPS_1_significant_gene_list.csv')
pps <- as.list(pps_df$'Genes')
pps_genes <- intersect(unlist(pps), seurat_genes)

taiji_df  <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/processed_Taiji/final_genes_all_sets_152.csv')
taiji <- as.list(taiji_df$'Genes')
taiji_genes <- intersect(unlist(taiji), seurat_genes)

###########################################################################

#Make a list of all gene lists of interest
###########################################################################
pathway_names <- c("PPI", "Taiji", "Review-2016", "Review-2022")
gene_set <- list(
  PPI = pps_genes,
  Taiji = taiji_genes,
  Review_2016 = vinuesa_genes,
  Review_2022 = tfh_review_genes
)
###########################################################################

#Run ssGSEA
###########################################################################
ssg_seurat = runEscape(seurat_obj, method = "ssGSEA", gene.sets = gene_set, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)

ssg_seurat = runEscape(subset_seurat, method = "ssGSEA", gene.sets = gene_set, min_size=2,
                       new.assay.name = "escape.ssGSEA", normalize = TRUE, BPPARAM = SnowParam(workers = 2), alpha=0.75)

DefaultAssay(ssg_seurat)  <- "escape.ssGSEA"


# Replace with your desired pathway name
pathway_name <- "PPI"

# Store ssGSEA score in metadata
ssg_seurat[[pathway_name]] <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")[pathway_name, ]


FeaturePlot(ssg_seurat, features = pathway_name, reduction = "spatial") +
  scale_fill_viridis_c(option = "plasma") +
  ggtitle(paste("ssGSEA score for", pathway_name))


FeaturePlot(ssg_seurat, features = "PPI", reduction = "spatial")

FeaturePlot(
  ssg_seurat,
  features = "PPI",  # or your ssGSEA pathway name
  assay = "escape.ssGSEA",
  reduction = "spatial",
  cols = c("blue", "white", "red")
)

row.names(ssg_seurat@assays$escape.ssGSEA)


subset_seurat <- subset(ssg_seurat, subset = cell_type %in% c("T_follicular_helper", "T_CD4"))

subset_seurat$PPI_score <- GetAssayData(ssg_seurat, assay = "escape.ssGSEA", slot = "data")["PPI", ]

FeaturePlot(
  subset_seurat,
  features = "PPI_score",
  reduction = "spatial",
  cols = c("blue", "red")
)

FeaturePlot(
  ssg_subset,
  features = "PPI_score",
  reduction = "spatial",
  cols = c("blue", "red")
)

plot1 <- FeaturePlot(
  ssg_seurat,
  features = "Review-2022",  # or your ssGSEA pathway name
  reduction = "spatial",
  cols = c("cyan2","salmon")
)
subset_seurat$cell_type <- factor(subset_seurat$cell_type, levels = c("T_follicular_helper", "T_CD4"))

plot2 <- DimPlot(subset_seurat, reduction = 'spatial', group.by = "cell_type", label = FALSE)

plot1

ggsave(filename ="/ix3/djishnu/Alisa/Tfh/scRNA/spatial/Tfh_vs_CD4T_xy_2022_reverse.pdf",
       plot = plot1,
       device = "pdf",
       width = 4,
       height = 4)



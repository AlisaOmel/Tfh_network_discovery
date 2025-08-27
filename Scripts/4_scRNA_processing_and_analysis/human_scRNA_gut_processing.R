library(Seurat)
library(ggplot2)
library(harmony)
library(Azimuth)

scRNA_data <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/GutCD4.RDS")

#################################################################################
#Remove outlier cells
#################################################################################
plot <- DimPlot(scRNA_data, reduction = 'umap.harmony', label = T, label.size =3, group.by='final_annotations')
cells.located <- CellSelector(plot = plot)

seurat_obj <- subset(scRNA_data, cells = cells.located)

saveRDS(seurat_obj, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/cluster_subset_GUT.rds")
#################################################################################


#################################################################################
#Recluster using harmony and Azimuth
#################################################################################

lafy <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/cluster_subset_GUT.rds")

mito.genes <- grep("^MT-", rownames(lafy), value = TRUE)
ribo.genes <- grep("^RPS|^RPL", rownames(lafy), value = TRUE)
hemo.genes <- grep("^HB[ABG]", rownames(lafy), value = TRUE)
sexgenes <- c("XIST", "TSIX", "DDX3Y", "KDM5D", "USP9Y", "UTY")  # XIST (X-inactivation), DDX3Y, KDM5D (Y-linked)
sexgenes <- intersect(sexgenes, rownames(lafy))
SFTP.genes <- grep("^SFT[ABCD]", rownames(lafy), value = TRUE)

lafy <- subset_GUT
DefaultAssay(lafy) <- "RNA"
lafy <-  Seurat::NormalizeData(lafy)
lafy <- FindVariableFeatures(lafy, selection.method = "vst", nfeatures = 5000)
lafy <- JoinLayers(lafy)
newvargenespostjoin <- VariableFeatures(object = lafy)
newvargenespostjoin<- setdiff(newvargenespostjoin, c(mito.genes, ribo.genes, hemo.genes, sexgenes, SFTP.genes)) #<- I use this to take out certain genes
lafy <- ScaleData(lafy, features = newvargenespostjoin)
lafy <- RunPCA(lafy, pc.genes = newvargenespostjoin, npcs = 100)

lafy <- RunHarmony(lafy, "Batch")
#Harmony converged after  iterations
harmony.embeddings <- Embeddings(lafy, reduction = "harmony")
ElbowPlot(lafy, ndims = 100, reduction = "harmony") #<-- this is if you want to see the ElbowPlot with harmony dimensions
#lafy <- RunUMAP(object = lafy, reduction = "harmony", dims = 1:20, verbose = FALSE, n.neighbors = 50)
lafy <- FindNeighbors(object = lafy, dims = 1:20, verbose = FALSE, k.param = 50, reduction = "harmony")
lafy <- FindClusters(object = lafy, verbose = FALSE, resolution = 0.7)

#UMAPPlot(lafy, label = T, label.size =3, group.by='Batch') 
#UMAPPlot(lafy, label = T, label.size =3, group.by='final_annotations')

lafy <- JoinLayers(lafy)
all_markers <- FindAllMarkers(lafy)

saveRDS(lafy, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/recluster_subset_GUT.rds")
#################################################################################



#################################################################################
#Relabel with Azimuth 
#################################################################################
#lafy <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/recluster_subset_GUT.rds")

mylabels <- RunAzimuth(lafy, reference = "tonsilref")

saveRDS(mylabels, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_Azimuth.rds")
                    


#################################################################################                    
#Remove Cycling T cells:
#################################################################################
gut_my_clustering <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/labels_GUT_Azimuth.rds")
plot_1 <- UMAPPlot(gut_my_clustering, reduction='umap', label = T,  group.by='predicted.celltype.l2') 
cells.located <- CellSelector(plot = plot_1)

gut_my_clustering_no_cycling_Tcells <- subset(gut_my_clustering, cells = cells.located)
saveRDS(gut_my_clustering_no_cycling_Tcells, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_Azimuth_nocyclingTcells.rds")



#################################################################################
#Identify Tfh
#################################################################################
gut_my_clustering_no_cycling_Tcells <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/labels_GUT_Azimuth_nocyclingTcells.rds")

plot_1 <- UMAPPlot(gut_my_clustering_no_cycling_Tcells, reduction='umap', label = T)# ,  group.by='predicted.celltype.l2') 

cluster_2_cells <- WhichCells(gut_my_clustering_no_cycling_Tcells, idents = "2")
Tfh_subset <- subset(gut_my_clustering_no_cycling_Tcells, subset = predicted.celltype.l2 %in% c( "GC-Tfh-SAP","Tfh T:B border"))
Tfh_subset <- subset(Tfh_subset,  cells = cluster_2_cells)
UMAPPlot(Tfh_subset, reduction='umap', label = T ,  group.by='predicted.celltype.l2') 
ncol(Tfh_subset)
saveRDS(Tfh_subset, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/.rds")



#################################################################################
#Identify Tfh non-Tfh
#################################################################################
exclude_groups <- c("GC-Tfh-SAP", "Tfh T:B border", "GC-Tfh-OX40", "Tfh-Mem", "Tfh-LZ-GC", 'Tregs','Eff-Tregs-IL32')
# Get all group labels
all_groups <- unique(gut_my_clustering_no_cycling_Tcells$predicted.celltype.l2)

# Select only the groups to keep
groups_to_keep <- setdiff(all_groups, exclude_groups)

# Use WhichCells to get those cells
cells_to_keep <- WhichCells(gut_my_clustering_no_cycling_Tcells, 
                            expression = predicted.celltype.l2 %in% groups_to_keep)

# Subset the object
cd4_subset <- subset(gut_my_clustering_no_cycling_Tcells, cells = cells_to_keep)
saveRDS(cd4_subset, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/non-Tfh_or_Treg_Cells.rds")

mylabels <- RunAzimuth(cd4_subset, reference = "pbmcref")

saveRDS(mylabels, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_AzimuthPBMC_nonTfh_noTfhorTreg.rds")
UMAPPlot(mylabels, reduction='umap', label = T ,  group.by='predicted.celltype.l2') 
#################################################################################



#################################################################################
#Identify CM Subset
#################################################################################
cd4_subset <- readRDS('/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_AzimuthPBMC_nonTfh_noTfhorTreg.rds')
cd4_subset_final <- subset(cd4_subset, subset =predicted.celltype.l2 %in% c( "CD4 TCM","CD4 TEM"))
plot_1 <-UMAPPlot(cd4_subset_final, reduction='umap', label = T ,  group.by='predicted.celltype.l2') 
cells.located <- CellSelector(plot = plot_1)
cd4_subset_final <- subset(cd4_subset_final, cells = cells.located)
saveRDS(cd4_subset_final, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_AzimuthPBMC_nonTfh_CMEM.rds")
#################################################################################



#################################################################################
#Visualize Marker genes
#################################################################################

nonTfh_filtered <- readRDS('"/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/non-Tfh_or_Treg_Cells.rds"')

genes_of_interest_CM <- c("CCR7", "SELL")
genes_present_CM <- genes_of_interest_CM[genes_of_interest_CM %in% rownames(nonTfh_filtered@assays$SCT@data)]
nonTfh_filtered$composite_score_cm<- colSums(nonTfh_filtered@assays$SCT@data[genes_present_CM, , drop=FALSE])

genes_of_interest_Th17 <- c("CCR6", "IL23R", "RORC")#, "IL17A")
genes_present_Th17 <- genes_of_interest_Th17[genes_of_interest_Th17 %in% rownames(nonTfh_filtered@assays$SCT@data)]

nonTfh_filtered$composite_score_Th17<- colSums(nonTfh_filtered@assays$SCT@data[genes_present_Th17, , drop=FALSE])

plot_1 <- UMAPPlot(nonTfh_filtered, reduction='harmony', label = T)#, group.by='predicted.celltype.l2', label.size =3)
plot_2 <- FeaturePlot(nonTfh_filtered, features = 'composite_score_cm', cols = c("lightgray", "red"), reduction= 'umap')
plot_1 + plot_2

combined_plot <- plot_1 + plot_2

ggsave(
  filename = "/ix3/djishnu/Alisa/Tfh/scRNA/modules_gs/CM_composite_umap.pdf",
  plot = combined_plot,
  width = 12, height = 6
)

#################################################################################


#################################################################################
#Th17 Cell Identification
#################################################################################
cluster_5_cells <- WhichCells(nonTfh_filtered, idents = c("5"))
Th17_cells <- subset(nonTfh_filtered, cells = cluster_5_cells)
plot_1 <- UMAPPlot(subset(Th17_cells, cells=cluster_5_cells ), reduction='umap', label = T, label.size =3, repel=TRUE, group.by= 'predicted.celltype.l2') #group.by="predicted.celltype.l2",

saveRDS(Th17_cells, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/Gut_Th17_only_Azimuth_cluster5.rds")


#################################################################################
#Get a subset of both non-Tfh and Tfh for the scGSEA analysis
#################################################################################
# Load original and subset objects
original <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/labels_GUT_Azimuth.rds")
tfh <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/Cluster2_GC_TB_Tfh_Cells.rds")
#non_tfh <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_AzimuthPBMC_nonTfh_CMEM.rds")
non_tfh <- readRDS("/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/Gut_scRNA/Gut_Th17_only_Azimuth_cluster5.rds")
# Get cell names from subsets
cells1 <- colnames(tfh)
cells2 <- colnames(non_tfh)

# Combine all cells of interest
all_cells <- union(cells1, cells2)

# Subset the original object
filtered <- subset(original, cells = all_cells)

# Create a new metadata column with the cell type
filtered$cell_type <- NA
filtered$cell_type[colnames(filtered) %in% cells1] <- "Tfh"
filtered$cell_type[colnames(filtered) %in% cells2] <- "non-Tfh"

# Check assignment
table(filtered$cell_type)

saveRDS(filtered, file = "/ix/djishnu/Alisa/Tfh/Ribeiro_Collaboration/sc_Data/labels_GUT_Th17_Tfh_combo.rds")


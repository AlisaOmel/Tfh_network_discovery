library(matrixStats)
library(readxl)
library(stringr)
library(msigdbr)
library(ggplot2)
library(dplyr)

#PPI network to generate the random networks from:
netwk <-read.delim("/ix/djishnu/Alisa/Tfh/Network_analysis/data/HomoSapiens_binary_co_complex_Feb2023_1_ppr_0.4.txt", header = F, sep=' ')
hop.distance <- rep(1, nrow(netwk))
df <- cbind.data.frame(netwk, hop.distance)
names(df) <- c("Node1", "Node2", "hop.distance")

#PPI + GRN Gene List
selectedGenes <- read.delim("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/pps_taiji_unique_genes_df.txt")

#The msigdb database to test:
Modules <- unique(selectedGenes$uGenes)                        

h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
GS_NAME <- unique(h_gene_sets$gs_name)

M <- matrix(0, nrow = length(GS_NAME), ncol = 6)
colnames(M) <- c('PathwayName', "Genes", "Node1", "Random1", "Node2", "Random2")

for(i in 1:length(GS_NAME)){
  
  gs <- GS_NAME[i]
  sub_gene_set <- subset.data.frame(h_gene_sets, gs_name == gs)
  TNF <- sub_gene_set$human_gene_symbol
  
  TNFSubset1 <- subset(df, ((df$Node1 %in% Modules) & (df$Node2 %in% TNF)))
  TNFSubset2 <- subset(df, ((df$Node2 %in% Modules) & (df$Node1 %in% TNF)))
  
  TNFSubset12 <- cbind.data.frame(c(TNFSubset1$Node1, TNFSubset2$Node2), c(TNFSubset1$Node2, TNFSubset2$Node1), c(TNFSubset1$hop.distance, TNFSubset2$hop.distance))
  names(TNFSubset12) <- c("Node1", "Node2", "hop.distance")
  
  test1 <- subset(TNFSubset12, hop.distance <= 1)
  
  M[i, 3] <- length(unique(test1$Node2))
  M[i, 5] <- length(unique(test1$Node1))
  M[i, 2] <- length(sub_gene_set$gene_symbol)
  M[i, 1] <- gs
  
  Pathways <- c()
  MODs <- c()
  
  j <- 1
  for (j in 1:1000){ 
    
    Random <- sample(setdiff(unique(c(df$Node1, df$Node2)), c(Modules,TNF)), size = length(Modules))
    
    randomSubset1 <- subset(df, ((df$Node1 %in% TNF) & (df$Node2 %in% Random)))
    randomSubset2 <- subset(df, ((df$Node2 %in% TNF) & (df$Node1 %in% Random)))
    
    randomSubset12 <- cbind.data.frame(c(randomSubset1$Node1, randomSubset2$Node2), c(randomSubset1$Node2, randomSubset2$Node1), c(randomSubset1$hop.distance, randomSubset2$hop.distance))
    names(randomSubset12) <- c("Node1", "Node2", "hop.distance")
    
    testRandom1 <- subset(randomSubset12, hop.distance <= 1)
    
    Pathways[j] <- length(unique(testRandom1$Node1))
    MODs[j] <- length(unique(testRandom1$Node2))
    
  }
  
  M[i, 4] <- mean(Pathways)
  M[i, 6] <- mean(MODs)
}

write.csv(M, "/ix/djishnu/Alisa/Tfh/Network_analysis/results/results_patternb_hallmark_1000.csv")

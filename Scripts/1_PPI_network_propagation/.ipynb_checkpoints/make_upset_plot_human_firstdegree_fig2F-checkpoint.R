library(ggplot2)
library(ComplexUpset)
library(pheatmap)
library(dplyr)

#Load all datasets
##################################################################################################
patternb_df  <- read.csv("/ix/djishnu/Alisa/Tfh/patternb_df.csv")
taiji_df  <- read.csv("/ix/djishnu/Alisa/Tfh/taiji_df.csv")
patternb_taiji_df  <- read.csv("/ix/djishnu/Alisa/Tfh/patternb_taiji_df.csv")
#logFC_RNA <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/log_FC_rna_set_unique_genes_df.csv')
#pps_noprop <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/PPS_noprop_set_184_unique_genes_df.csv')

logFC_RNA <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/log_FC_rna_set_correct_184_0.01_unique_genes_df.csv')
pps_noprop <- read.csv('/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/PPS_noprop_set_184_correct_unique_genes_df.csv')

patternb <- as.list(patternb_df$'Genes')
taiji <- as.list(taiji_df$'Genes')
patternb_taiji <- as.list(patternb_taiji_df$'Genes')
logFC_RNA <- as.list(logFC_RNA$"Genes")
pps_noprop <- as.list(pps_noprop$"Genes")


#Randoms were made with a python code  ________ for 1000 permutations 
#patternb_random_df <- read.csv("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/pps_fd_random_upset_human.csv")
#patternb_random <- as.list(patternb_random_df$'Genes')
patternb_random_df <- read.csv("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/pps_randomfd_random_upset_human.csv")
patternb_random <- as.list(patternb_random_df$'Genes')

#taiji_random_df <- read.csv("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/taiji_fd_random_upset_human.csv")
#taiji_random <- as.list(taiji_random_df$'Genes')

taiji_random_df <- read.csv("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/taiji_randomfd_random_upset_human.csv")
taiji_random <- as.list(taiji_random_df$'Genes')

#patternb_taiji_random_df <- read.csv("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/pps_taiji_fd_random_upset_human.csv")
#patternb_taiji_random <- as.list(patternb_taiji_random_df$'Genes')

patternb_taiji_random_df <- read.csv("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/pps_taiji_randomfd_random_upset_human.csv")
patternb_taiji_random <- as.list(patternb_taiji_random_df$'Genes')

logFC_random_df <- read.csv("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/logFC_noprop_randomfd_random_upset_human.csv")
logFC_random <- as.list(logFC_random_df$'Genes')

pps_noprop_random_df <- read.csv("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/pps_noprop_randomfd_random_upset_human.csv")
pps_noprop_random <- as.list(pps_noprop_random_df$'Genes')

#logFC_random_df <- read.csv("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/logFC_RNA_fd_correct_184_0.01_random_upset_human.csv")
#logFC_random <- as.list(logFC_random_df$'Genes')

#pps_noprop_random_df <- read.csv("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/PPS_noprop_fd_184_correct_random_upset_human.csv")
#pps_noprop_random <- as.list(pps_noprop_random_df$'Genes')



vinuesa_df  <-read.csv("/ix/djishnu/Alisa/Tfh/vinuesa_list.csv")
tfh_review_df  <- read.csv("/ix/djishnu/Alisa/Tfh/TFH_review_list.csv")
vinuesa <- as.list(vinuesa_df$'Genes')
tfh_review <- as.list(tfh_review_df$'Genes')
##################################################################################################

#Choose the datasets to be combined, the upset plot needs this format: 
##################################################################################################
combine_gene_lists <- function(gene_lists, column_names) {
  # Ensure the number of column names matches the number of gene lists
  if (length(gene_lists) != length(column_names)) {
    stop("The number of column names must match the number of gene lists.")
  }
  
  # Get unique genes from all lists
  unique_genes <- unique(unlist(gene_lists))
  
  # Initialize a data frame with the unique genes
  df <- data.frame(Genes = unique_genes)
  
  # Add columns for each gene list indicating presence (1) or absence (0)
  for (i in seq_along(gene_lists)) {
    df[[column_names[i]]] <- as.integer(df$Genes %in% gene_lists[[i]])
  }
  
  return(df)
}

#Network Propagation graph
gene_lists<- list(patternb_taiji_random, taiji_random, patternb_random, vinuesa, tfh_review, patternb,  taiji, patternb_taiji)
#colname <- c("PPI+GR Random", "GR Random", "PPI Random", "Vinuesa et al.", "Hart et al.", "PPI", "Taiji", "PPS+Taiji", 'logFC_RNA')
colname <- c("PPI+GRN Random", "GRN Random", "PPI Random", "2016 Review", "2022 Review", "PPI", "GRN", "PPI+GRN")

#gene_lists <- c(taiji_random, taiji)
#colname <- c("GRN Random", "GRN", "2016 Review", "2022 Review")
#Just scores no network prop graph
#gene_lists<- list(pps_noprop_random, logFC_random, vinuesa, tfh_review, logFC_RNA, pps_noprop)
#colname <- c("PPI+GR Random", "GR Random", "PPI Random", "Vinuesa et al.", "Hart et al.", "PPI", "Taiji", "PPS+Taiji", 'logFC_RNA')
#colname <- c("PPS Random", "logFC RNA Random", "2016 Review", "2022 Review", 'logFC RNA', "PPS")
##################################################################################################


#Calculate P-values
##################################################################################################
result <- combine_gene_lists(gene_lists, colname)
df <- subset(result, select = -Genes)

random_sizes <- list(
  "PPI Random" = 3626,
  "GRN Random" = 3147,
  "PPI+GRN Random" = 5350
)

#random_sizes <- list(
#  "PPS Random" = 3626,
#  "logFC RNA Random" = 3626
#)


#compare_to_random <- function(df, method_col, random_col, reference_col) {
#  a1 <- sum(df[[method_col]] == 1 & df[[reference_col]] == 1)
#  m  <- sum(df[[method_col]] == 1)
  
#  a2 <- sum(df[[random_col]] == 1 & df[[reference_col]] == 1)
#  r  <- sum(df[[random_col]] == 1)
  
#  mat <- matrix(c(a1, m - a1, a2, r - a2), nrow = 2, byrow = TRUE)
#  fisher.test(mat, alternative = "greater")$p.value
#}

#compare_to_random <- function(df, method_col, random_col, reference_col) {
#  # Overlaps with reference
#  a1 <- sum(df[[method_col]] == 1 & df[[reference_col]] == 1)
#  m  <- sum(df[[method_col]] == 1)
  
#  a2 <- sum(df[[random_col]] == 1 & df[[reference_col]] == 1)
#  r  <- 3147#sum(df[[random_col]] == 1)

  
#  # Perform binomial proportions test
#  test <- prop.test(x = c(a1, a2), n = c(m, r), alternative = "greater", correct = FALSE)
#  return(test$p.value)
#}

compare_to_random <- function(df, method_col, random_col, reference_col, random_sizes) {
  a1 <- sum(df[[method_col]] == 1 & df[[reference_col]] == 1)
  m  <- sum(df[[method_col]] == 1)
  print(m)
  
  a2 <- sum(df[[random_col]] == 1 & df[[reference_col]] == 1)
  r  <- random_sizes[[random_col]]  # Use average value from list
  
  test <- prop.test(x = c(a1, a2), n = c(m, r), alternative = "greater", correct = FALSE)
  return(test$p.value)
}

# --- Define your columns ---
methods   <- c("GRN", "PPI", "PPI+GRN")
#methods   <- c( "PPS", 'logFC RNA')
randoms   <- c( "GRN Random", "PPI Random", "PPI+GRN Random")
#randoms   <- c("PPS Random", "logFC RNA Random")
references <- c("2016 Review", "2022 Review")
# --- Initialize p-value matrix ---
pval_matrix <- matrix(nrow = length(methods), ncol = length(references),
                      dimnames = list(methods, references))

# --- Calculate p-values ---
#for (i in seq_along(methods)) {
#  print(i)
#  method_col <- methods[i]
#  random_col <- randoms[i]
  
#  for (ref in references) {
#    pval_matrix[method_col, ref] <- compare_to_random(df, method_col, random_col, ref)
#  }
#}

for (i in seq_along(methods)) {
  method_col <- methods[i]
  random_col <- randoms[i]
  
  for (ref in references) {
    pval_matrix[method_col, ref] <- compare_to_random(df, method_col, random_col, ref, random_sizes)
  }
}



# --- Save to CSV ---
write.csv(pval_matrix, file = "/ix/djishnu/Alisa/Tfh/method_vs_random_pvalues_figF_network.csv", quote = FALSE)
##################################################################################################

#Plot network propagation upset plot
##################################################################################################
p <- upset(
  df, colname, 
  mode = 'inclusive_intersection',
  min_size = 0,
  intersections = list(
    c('PPI+GRN', '2022 Review'),
    c('PPI+GRN Random', '2022 Review'),
    c('GRN', '2022 Review'),
    c('GRN Random', '2022 Review'),
    c('PPI', '2022 Review'),
    c('PPI Random', '2022 Review'),
    c('PPI+GRN', '2016 Review'),
    c('PPI+GRN Random', '2016 Review'),
    c('GRN', '2016 Review'),
    c('GRN Random', '2016 Review'),
    c('PPI', '2016 Review'),
    c('PPI Random', '2016 Review')
  ),
  #sort_sets = FALSE,
  sort_intersections = FALSE,
  base_annotations=list('Intersection size'=intersection_size( mode='inclusive_intersection')), #counts=FALSE, 
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(text = element_text(size = 20)),
      'overall_sizes' = theme(text = element_text(size = 20)),
      'set_sizes' = theme(text = element_text(size = 20)),
      'set_names' = theme(text = element_text(size = 20)),
      'intersection_size' = theme(  # ← This is key for the y-axis
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18)
      ),
      'default' = theme(text = element_text(size = 20))
    )
  ),
  
  set_size = FALSE,
  queries = list(
    upset_query(intersect = c('PPI+GRN', '2016 Review'), color = "olivedrab3", fill = 'olivedrab3'),
    upset_query(intersect = c('PPI+GRN', '2022 Review'), color = "olivedrab3", fill = "olivedrab3"),
    upset_query(intersect = c('PPI+GRN Random', '2016 Review'), color = "grey80", fill = 'grey80'),
    upset_query(intersect = c('PPI+GRN Random', '2022 Review'), color = "grey80", fill = "grey80"),
    upset_query(intersect = c('PPI', '2016 Review'), color = "steelblue2", fill = "steelblue2"),
    upset_query(intersect = c('PPI', '2022 Review'), color = "steelblue2", fill = "steelblue2"),
    upset_query(intersect = c('GRN', '2022 Review'), color = "gold", fill = 'gold'),
    upset_query(intersect = c('GRN', '2016 Review'), color = "gold", fill = "gold"),
    upset_query(intersect = c('PPI Random', '2016 Review'), color = "grey80", fill = "grey80"),
    upset_query(intersect = c('PPI Random', '2022 Review'), color = "grey80", fill = "grey80"),
    upset_query(intersect = c('GRN Random', '2022 Review'), color = "grey80", fill = 'grey80'),
    upset_query(intersect = c('GRN Random', '2016 Review'), color = "grey80", fill = "grey80")
  )
)

p



#pdf("/ix/djishnu/Alisa/Tfh/ForPaper/upset_plots/human_fd_final_correctrandom.pdf", width = 8.5, height = 10)
pdf("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/human_fd_final_correctrandom.pdf", width = 8.5, height = 10)
print(p)
dev.off()

# Save as PNG
#png("/ix/djishnu/Alisa/Tfh/ForPaper/upset_plots/human_fd_final_correctrandom.png", width = 8.5, height = 10, units = "in", res = 300)
png("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/human_fd_final_correctrandom.png", width = 8.5, height = 10, units = "in", res = 300)
print(p)
dev.off()
##################################################################################################

#Plot just scores plot
##################################################################################################
p <- upset(
  df, colname, 
  mode = 'inclusive_intersection',
  min_size = 0,
  intersections = list(
    c('logFC RNA', '2022 Review'),
    c('logFC RNA Random', '2022 Review'),
    c("PPS", "2022 Review"),
    c('PPS Random', "2022 Review"),
    c('logFC RNA', '2016 Review'),
    c('logFC RNA Random', '2016 Review'),
    c('PPS', "2016 Review"),
    c('PPS Random', "2016 Review")
  ),
  #sort_sets = FALSE,
  sort_intersections = FALSE,
  base_annotations=list('Intersection size'=intersection_size( mode='inclusive_intersection')), #counts=FALSE, 
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(text = element_text(size = 20)),
      'overall_sizes' = theme(text = element_text(size = 20)),
      'set_sizes' = theme(text = element_text(size = 20)),
      'set_names' = theme(text = element_text(size = 20)),
      'intersection_size' = theme(  # ← This is key for the y-axis
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18)
      ),
      'default' = theme(text = element_text(size = 20))
    )
  ),
  
  set_size = FALSE,
  queries = list(
    upset_query(intersect = c('PPS', '2016 Review'), color = "#DDA0DD" , fill = "#DDA0DD"),
    upset_query(intersect = c('PPS', '2022 Review'), color = "#DDA0DD", fill = "#DDA0DD"),
    upset_query(intersect = c('PPS Random', '2016 Review'), color = "grey80", fill = "grey80"),
    upset_query(intersect = c('PPS Random', '2022 Review'), color = "grey80", fill = "grey80"),
    upset_query(intersect = c('logFC RNA', '2016 Review'), color = "#F08080", fill = "#F08080"),
    upset_query(intersect = c('logFC RNA', '2022 Review'), color = "#F08080", fill = "#F08080"),
    upset_query(intersect = c('logFC RNA Random', '2022 Review'), color = "grey80", fill = 'grey80'),
    upset_query(intersect = c('logFC RNA Random', '2016 Review'), color = "grey80", fill = "grey80")
    
  )
)

p

pdf("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/human_only_w_noprop_controls_correct_184_0.01.pdf", width = 5.67, height = 10)
print(p)
dev.off()

# Save as PNG
png("/ix3/djishnu/Alisa/Tfh/upsetplot_randoms/human_fd_only_noprop_controls_correct_184_0.01.png", width = 5.67, height = 10, units = "in", res = 300)
print(p)
dev.off()
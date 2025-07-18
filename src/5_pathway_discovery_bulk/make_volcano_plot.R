library(matrixStats)
library(readxl)
library(ggplot2)
library(dplyr)

data <-read.delim("/ix3/djishnu/Alisa/Tfh/mouse_bulk/results_patternb_mouse_BIOCARTA_1000.csv", sep=',')


#######################################################################################
#Q-value calculation
#######################################################################################

pVAL <- c()
for(i in 1:nrow(data)){
  # test <- prop.test(x = c(data$V1[i], data$V2[i]), n = c(data$Genes[i], data$Genes[i]), alternative = "greater")
  test <- prop.test(x = c(data$Node1[i], data$Random1[i]), n = c(data$Genes[i], data$Genes[i]), correct= FALSE)
  pVAL[i] <- test$p.value
}

adj_pval = p.adjust(pVAL, method = "BH")
print(adj_pval)

logQ <- -1*log10(adj_pval)

#######################################################################################
#Calculate fold change using odds ratio, secondary version adding epsilon to reduce Inf
#######################################################################################

df <- data[, c(2,3,4,5,6)]
df <- cbind.data.frame(df, pVAL)
f1 <- (df$Node1)/df$Genes
f2 <- (df$Random1)/df$Genes
#oddsRatio <- (f1 /(1- f1)) / (f2 /(1 - f2))
#fc <- log2(oddsRatio)

epsilon <- 1e-5
oddsRatio <- ((f1 + epsilon) / (1 - f1 + epsilon)) / ((f2 + epsilon) / (1 - f2 + epsilon))
fc <- log2(oddsRatio)

#######################################################################################

#######################################################################################
#Make dataframes for plot
#######################################################################################

pathwayName <- df[, c(1)]

#plotDF <- cbind.data.frame(pathwayName, fc, P, logP)
plotDF <- cbind.data.frame(pathwayName, fc, adj_pval, logQ)

plotDF <- na.omit(plotDF)
volcano = ggplot(data = plotDF, aes(x = fc, y = logQ))
volcano + geom_point() + xlab(expression(log[2](FC))) + ylab(expression(log[10](Q))) + xlim(-2,2) + ggtitle("pathway vs random") +
  theme(plot.title = element_text(hjust = 0.5))


myDiff1p <- plotDF

myDiff1p <- myDiff1p %>%
  mutate(threshold = factor(case_when(fc > 0.58 & logQ > 0.7 ~ "cond1",
                                      fc < -0.58 & logQ > 0.7 ~ "cond2",
                                      TRUE ~ "cond3")))
write.csv(myDiff1p, "/ix3/djishnu/Alisa/Tfh/mouse_bulk/patternb_taiji_BIOCARTA_1000_volc_plot_thresh_qval.csv")

#######################################################################################



#######################################################################################
#Plotting
#######################################################################################

volcano <- ggplot(data=myDiff1p, aes(x = fc, y = logQ)) + 
  geom_point(aes(color = threshold), alpha=1, size = 1.75) + 
  geom_vline(xintercept=c(-0.58, 0.58), color="olivedrab3", alpha=0.3) +
  geom_hline(yintercept= 0.70, color="steelblue2", alpha=0.5) +
  xlab(expression(log[2](odds_ratio))) + 
  ylab(expression(-log[10](Q))) + 
  theme_classic() +
  xlim(c(-3, 3)) +
  ylim(c(0, 3)) +
  scale_color_manual(name = "Threshold",
                     values = c("cond1" = "steelblue2", "cond2" = "black", "cond3" = "black")) +
  theme(
    axis.text = element_text(size = 14),   # Adjust size of axis numbers
    axis.title = element_text(size = 16), # Adjust size of axis labels
    legend.position = "None"
  ) +
  annotate("text", x = c(2, 2, 2, 2), y = c(2.95, 1.8, 1.25, 0.9), size = 4.5, label = c('BIOCARTA_MAPK_PATHWAY',"BIOCARTA_NFAT_PATHWAY", "BIOCARTA_IL6_PATHWAY", "BIOCARTA_IL12_PATHWAY"), color = "steelblue2", fontface = 2)

volcano

ggsave(
  filename = "/ix3/djishnu/Alisa/Tfh/mouse_bulk/BIOCARTA_volcano_plot.pdf",
  plot = volcano,
  width = 6,
  height = 6 
)

dev.off()
#######################################################################################
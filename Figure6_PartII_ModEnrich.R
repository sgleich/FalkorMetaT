### Figure 6; Part II- Module enrichment analyses ###
### This script will allow you to perform module enrichment analyses using the clusterProfiler package. This analysis uses the results obtained a differential expression analysis (such as that carried out in Figure6_PartI_DinoDE.R) ###
### Last Updated: July 31, 2023 ###

# Load libraries
library(clusterProfiler)
library(stringr)
library(tidyverse)

# Load DGE analysis output
dge <- read.csv("DE_150m_Falkor.csv",header=TRUE)
dge$KO_list <- str_remove_all(dge$KO_list,"ko:")
dge$KO_list <- str_split(dge$KO_list,"_")
dge<- unnest(dge,KO_list)

# Subset for orthologous groups that were 1.) upregulated in the cyclone relative to the anticyclone and had an adjusted p-value less than 0.01 AND 2.) that were upregulated in the anticyclone relative to the cyclone and had an adjusted p-value less than 0.01.
cy <- subset(dge, logFC > 0 & padj < 0.01)
acy <- subset(dge, logFC < 0 & padj < 0.01)

# Background gene list
bg <- unique(dge$KO_list)

# Run enrichMKEGG function on cy and acy dataframes
cyOut <- enrichMKEGG(cy$KO_list, pvalueCutoff=0.05,organism="ko",pAdjustMethod="fdr",keyType='kegg',universe=bg)
cyOut <- cyOut@result

acyOut <- enrichMKEGG(acy$KO_list, pvalueCutoff=0.05,organism="ko",pAdjustMethod="fdr",keyType='kegg',universe=bg)
acyOut <- acyOut@result

# Only plot modules that had a clusterProfiler adjusted p-value < 0.05
finCy <- subset(cyOut,p.adjust<0.05)
finAcy <- subset(acyOut,p.adjust<0.05)

# Add eddy info and combine dataframes
finCy$Eddy <- "Cyclonic"
finAcy$Eddy <- "Anticyclonic"
full <- rbind(finCy,finAcy)

# Plot results
p <- ggplot(full,aes(x=Count,y=reorder(Description,-Count),fill=Eddy))+geom_bar(stat="identity",width = 0.9, position = position_dodge(),color="black")+theme_classic(base_size = 10)+theme(legend.position = "right")+ylab("KEGG Module")+ggtitle("Dinoflagellate Module Enrichment at 150 m")+geom_vline(xintercept=0,linetype="dashed")+scale_fill_manual(values=c("red","blue"),breaks=c("Anticyclonic","Cyclonic"))+xlab("Number of Hits to KEGG Module")+xlim(0,21)
# ggsave("FigureX.pdf",width=15,height=8)

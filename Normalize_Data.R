### This script will allow you to normalize transcript abundance data across samples using the edgeR TMM CPM method. ###
### By: Samantha Gleich ###
### Last Updated: November 10, 2023 ###

# Load libraries
library(edgeR)

# Load in wide format data (compiled but not yet normalized; obtained via Compile_Data.R)
df <- read.csv("wide_format_data.csv",header=TRUE,row.names=1)
df <- subset(df,grepl("Eukaryota",df$Taxonomy))
df$KEGG <- as.character(df$KEGG)

df$KEGG <- ifelse(is.na(df$KEGG),"",df$KEGG)

# Set up edgeR list. Counts are salmon sample counts, genes are identifiers (i.e. Name, KEGG, and Taxonomy), groups are you telling edgeR which groups exist in your dataset (i.e. sample types)
dge.metaT <- DGEList(counts=df[4:27],genes=df[1:3],group=c(rep("Cyclonic-DCM",3),rep("Cyclonic-25m",3),rep("Anticyclonic-250m",3),rep("Anticyclonic-150m",3),rep("Anticyclonic-DCM",3),rep("Anticyclonic-25m",3),rep("Cyclonic-250m",3),rep("Cyclonic-150m",3)))

# Optional check to see if groupings look good
# dge.metaT$samples

# Use TMM normalization method in edgeR
metaT.norm <- calcNormFactors(dge.metaT,method = "TMM")
cpm_data <- cpm(metaT.norm, normalized.lib.sizes=TRUE, log=FALSE) 
cpm_data<-as.data.frame(cpm_data)  
data_CPM<-data.frame(metaT.norm$genes,cpm_data)

# Save normalized dataset
write.csv(data_CPM, "normalized_metaT_data.csv")

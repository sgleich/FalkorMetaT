### Figure 3 - NMDS/PCA ###
### This script will allow you to make an NMDS plot using rRNA reads pull out of metatranscriptomes and a PCA plot using transcripts from metatranscriptomes. ###
### Last Updated: August 1, 2023 ###

# Load libraries 
library(ggplot)
library(vegan)
library(tidyverse)
library(reshape2)
library(patchwork)

### NMDS PLOT ###
# Compile miTag results (taxonomy associated with rRNA reads from sortmerna)
samples <- c("10-cyclonic-112","11-cyclonic-112","12-cyclonic-112","13-cyclonic-25","14-cyclonic-25","15-cyclonic-25","19-anticyclonic-250","20-anticyclonic-250","21-anticyclonic-250","22-anticyclonic-150","23-anticyclonic-150","24-anticyclonic-150","25-anticyclonic-122","26-anticyclonic-122","27-anticyclonic-122","28-anticyclonic-25","29-anticyclonic-25","30-anticyclonic-25","4-cyclonic-250","5-cyclonic-250","6-cyclonic-250","7-cyclonic-150","8-cyclonic-150","9-cyclonic-150")

# Loop to read in all samples
all <- NULL
for (sample in samples){
  df <- read.delim(paste(sample,"tax_assignments.txt",sep="_"),header=FALSE,row.names=NULL)
  colnames(df)<- c("seqID","Tax","X1","X2")
  df$Sample <- sample
  df$Count <- 1
  all <- rbind(all,df)}

# Remove columns we don't need
all$seqID <- NULL
all$Num1 <- NULL
all$Num2 <- NULL

# Find total number of reads associated with each taxonomic ID and each sample
allSum <- all %>% group_by(Tax,Sample) %>% summarize(sum(Count)) %>% as.data.frame()

# Remove unassigned miTags - we will only consider seqs that were annotated at or below "Eukaryote"
mitag <- subset(allSum,Tax!="Unassigned")

# Turn long dataframe into wide dataframe
mitag <- mitag %>% pivot_wider(id_cols=Sample,names_from = Tax,values_from = sum.Count.,values_fill = 0) %>% as.data.frame()
rownames(mitag) <- mitag$Sample
mitag$Sample <- NULL

# Carry out NMDS
mitagNorm <- decostand(mitag,method="total")
rowSums(mitagNorm) # Row sums should equal 1

mitagBray <- vegdist(mitagNorm,method="bray") # Bray-curtis dissimiliarity 
nmdsDf <- metaMDS(mitagBray,k=2, distance="bray") # NMDS in 2 dimensions
nmdsDf$stress # Record stress: 0.059

# Extract results of NMDS
nmdsDf <- nmdsDf$points
cols <- colsplit(rownames(nmdsDf),"-",c("num","eddy","depth"))
nmdsDf <- as.data.frame(nmdsDf)
nmdsDf$eddy <- cols$eddy
nmdsDf$depth <- as.factor(cols$depth)

# Plot NMDS
nmdsPlot <- ggplot(nmdsDf,aes(MDS1,MDS2,shape=eddy))+geom_point(size=3,aes(fill=depth))+scale_shape_manual(values=c(22,21),name="Eddy",labels=c("Anticyclonic","Cyclonic"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("NMDS1")+ylab("NMDS2")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("Community composition (18S rRNA)")+scale_fill_manual(values=c("darkolivegreen3","dodgerblue3","indianred","darkgoldenrod1"),name="Depth",labels=c("25 m","DCM","150 m","250 m"))

### PCA PLOT ###
# Load in normalized dataframe (obtained via Normalize_Data.R)
df <- read.csv("normalized_metaT_data.csv",header=TRUE,row.names=1)

# Remove columns we don't need
df$contigID <- NULL
df$Taxonomy <- NULL
df$Cluster <- NULL 

# Remove rows without KO term
df <- subset(df,KO!="")

# Sum counts per KO terms for each sample
dfSum <- df %>% group_by(KO) %>% summarize_all(sum) %>% as.data.frame()

# Run PCA
rownames(dfSum) <- dfSum$KO
dfSum$KO <- NULL
dfSum <- as.data.frame(t(dfSum))
pcaOut <- prcomp(dfSum,scale=TRUE)

# Extract results of PCA
pcaPlot <- data.frame(pcaOut$x)
colz <- colsplit(rownames(pcaPlot),"\\.",c("x","eddy","depth","rep"))
pcaPlot$eddy <- colz$eddy
pcaPlot$depth <- colz$depth

# See percent of variance explained by the first 2 axes
summary(pcaOut)

# Plot PCA
pcaPlot <- ggplot(pcaPlot,aes(PC1,PC2,shape=eddy))+geom_point(size=3,aes(fill=depth))+scale_shape_manual(values=c(22,21),name="Eddy",labels=c("Anticyclonic","Cyclonic"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("PC1 (29.7%)")+ylab("PC2 (16.6%)")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("Metabolic potential (mRNA)")+scale_fill_manual(breaks=c("25m","112m","150m","250m"),values=c("darkolivegreen3","dodgerblue3","indianred","darkgoldenrod1"),name="Depth",labels=c("25 m","DCM","150 m","250 m"))

# Combine plots with a common legend and add panel labels (a-b)
nmdsPlot+ pcaPlot + plot_layout(guides = "collect",nrow=2)+plot_annotation(tag_levels="a")
# ggsave("Figure3.pdf",width=6,height=8)

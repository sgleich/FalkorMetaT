### Figure 3 - NMDS/PCA ###
### This script will allow you to make an NMDS plot using rRNA reads pull out of metatranscriptomes and a PCA plot using transcripts from metatranscriptomes. ###
### Last Updated: July 23, 2023 ###

# Load libraries 
library(ggplot)
library(vegan)
library(tidyverse)
library(reshape2)
library(patchwork)

# NMDS Plot
mitag <- read.csv("miTag_final_april2023.csv",header=TRUE,row.names=1)
mitag <- subset(mitag,Tax!="Unassigned")

mitag <- mitag %>% pivot_wider(id_cols=Sample,names_from = Tax,values_from = sum.Count.,values_fill = 0) %>% as.data.frame()
rownames(mitag) <- mitag$Sample
mitag$Sample <- NULL

mitagNorm <- decostand(mitag,method="total")
rowSums(mitagNorm)
mitagBray <- vegdist(mitagNorm,method="bray") # Bray-curtis dissimiliarity 
nmdsDf <- metaMDS(mitagBray,k=2)
nmdsDf$stress # Record stress - 0.059
nmdsDf <- nmdsDf$points
colz <- colsplit(rownames(nmdsDf),"-",c("num","eddy","depth"))
nmdsDf <- as.data.frame(nmdsDf)
nmdsDf$eddy <- colz$eddy
nmdsDf$depth <- as.factor(colz$depth)

nmdsPlot <- ggplot(nmdsDf,aes(MDS1,MDS2,shape=eddy))+geom_point(size=3,aes(fill=depth))+scale_shape_manual(values=c(22,21),name="Eddy",labels=c("Anticyclonic","Cyclonic"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("NMDS1")+ylab("NMDS2")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("Community composition (18S rRNA)")+scale_fill_manual(values=c("darkolivegreen3","dodgerblue3","indianred","darkgoldenrod1"),name="Depth",labels=c("25 m","DCM","150 m","250 m"))
nmdsPlot

# PCA Plot
pcaDf <- read.csv("KOSums_AllEuk_May2023.csv",header=TRUE,row.names=1)
rownames(pcaDf) <- pcaDf$KO
pcaDf$KO <- NULL
pcaDf <- as.data.frame(t(pcaDf))
pcaOut <- prcomp(pcaDf,scale=TRUE)
pcaPlot <- data.frame(pcaOut$x)
colz <- colsplit(rownames(pcaPlot),"\\.",c("x","eddy","depth","rep"))
pcaPlot$eddy <- colz$eddy
pcaPlot$depth <- colz$depth

summary(pcaOut)

pcaP <- ggplot(pcaPlot,aes(PC1,PC2,shape=eddy))+geom_point(size=3,aes(fill=depth))+scale_shape_manual(values=c(22,21),name="Eddy",labels=c("Anticyclonic","Cyclonic"))+guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+xlab("PC1 (29.7%)")+ylab("PC2 (16.6%)")+geom_vline(xintercept = 0,linetype="dotted")+geom_hline(yintercept = 0,linetype="dotted")+ggtitle("Metabolic potential (mRNA)")+scale_fill_manual(breaks=c("25m","112m","150m","250m"),values=c("darkolivegreen3","dodgerblue3","indianred","darkgoldenrod1"),name="Depth",labels=c("25 m","DCM","150 m","250 m"))
pcaP

nmdsPlot+ pcaP + plot_layout(guides = "collect",nrow=2)+plot_annotation(tag_levels="a")
# ggsave("Figure3.pdf",width=6,height=8)

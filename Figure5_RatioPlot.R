### Figure 5 - Phototrophy:Heterotrophy ratio plots ###
### This script will allow you to calculate and plot phototrophy:heterotrophy ratios from a list of phototrophy-associated and heterotrophy-associated KO terms that are pre-defined ###
### Last Updated: November 10, 2023 ###

# Load libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(patchwork)

# Load in normalized dataframe (obtained via Normalize_Data.R)
# Here we will load in the normalized dataframe that contains all contigs annotated as "Eukaryota". To run this script for other groups, you will want to normalize the dataframe for specific taxa (e.g., dinoflagellates, diatoms, etc.).
all <- read.csv("normalized_metaT_data.csv",header=TRUE,row.names=1)

# Remove columns we don't need
all$Name <- NULL
all$Taxonomy <- NULL

# Remove rows without KO term
all <- subset(all,KEGG!="")

# Sum counts per KO terms for each sample
dfSum <- all %>% group_by(KEGG) %>% summarize_all(sum) %>% as.data.frame()

# Ratio plot function
# Input: dataframe and plot title
# Output: Plot showing phototrophy:heterotrophy transcript ratios +/- SE across depths and eddy types
makeRatioPlot <- function(df,title){
  df <- subset(df,!is.na(KEGG))
  df <- subset(df,KEGG !="")
  
  # Gene groups
  `%ni%` <- Negate(`%in%`)
  ko <- read.csv("KEGG.csv",header=TRUE)
  
  photo <- c("Photosynthesis","Calvin cycle","Antenna proteins")
  
  koPhoto <- subset(ko,Biomarker.Gene.Group %in% photo)  
  koHetero <- subset(ko,Biomarker.Gene.Group %ni% photo) 
  koHetero <- subset(koHetero,Biomarker.Gene.Group!="Inorganic N uptake and assimilation")
  # print(unique(koHetero$Biomarker.Gene.Group))
  
  dfPhoto <- subset(df,KO %in% koPhoto$KEGG.KO)
  dfHetero <- subset(df,KO %in% koHetero$KEGG.KO)
  
  dfPhotoMelt <- melt(dfPhoto,id.vars="KO")
  dfHeteroMelt <- melt(dfHetero,id.vars="KO")
  
  colzPhoto <- colsplit(dfPhotoMelt$variable,"\\.",c("Num","Eddy","Depth","Rep"))
  colzHetero <- colsplit(dfHeteroMelt$variable,"\\.",c("Num","Eddy","Depth","Rep"))
  
  dfPhotoMelt$Depth <- colzPhoto$Depth
  dfHeteroMelt$Depth <- colzHetero$Depth
  dfPhotoMelt$Eddy <- colzPhoto$Eddy
  dfHeteroMelt$Eddy <- colzHetero$Eddy
  dfPhotoMelt$Rep <- colzPhoto$Rep
  dfHeteroMelt$Rep <- colzHetero$Rep
  
  dfPhotoMelt <- dfPhotoMelt %>% group_by(Depth,Eddy,Rep) %>% summarize(sPhoto=sum(value)) %>% as.data.frame()
  dfHeteroMelt <- dfHeteroMelt %>% group_by(Depth,Eddy,Rep) %>% summarize(sHetero=sum(value)) %>% as.data.frame()
  
  dfFull <- left_join(dfPhotoMelt,dfHeteroMelt)
  dfFull$Ratio <- dfFull$sPhoto/dfFull$sHetero
  
  dfFin <- dfFull%>%group_by(Eddy,Depth)%>%summarize(meanz=mean(Ratio),sez=(sd(Ratio)/sqrt(3)))
  
  dfFin$Depth <- ifelse(dfFin$Depth=="112m","DCM", dfFin$Depth)
  dfFin$Depth <- ifelse(dfFin$Depth=="25m","25 m", dfFin$Depth)
  dfFin$Depth <- ifelse(dfFin$Depth=="150m","150 m", dfFin$Depth) 
  dfFin$Depth <- ifelse(dfFin$Depth=="250m","250 m", dfFin$Depth)
  dfFin$Depth <- factor(dfFin$Depth,levels=c("250 m","150 m","DCM","25 m"))
  
  p <- ggplot(dfFin, aes(x=Depth,y=meanz,fill=Eddy,group=Eddy))+geom_bar(stat='identity', position='dodge',color="black")+scale_y_continuous(breaks = scales::pretty_breaks(n = 3),position='right')+geom_errorbar(aes(ymin=meanz-sez, ymax=meanz+sez), width=0.4,position=position_dodge(.9))+coord_flip()+theme_classic(base_size = 14)+xlab("Depth (m)")+ggtitle(paste(title))+scale_fill_manual(values=c("red","blue"))+ylab("Phototrophy:Heterotrophy Transcript Ratio +/- SE")
  return(p)}

# Run makeRatioPlot function for all different taxa
allEuk <- makeRatioPlot(dfSum,"All eukaryotes")
# chlorophyte <- makeRatioPlot(chloro,"Chlorophyte")
# ciliate <- makeRatioPlot(cil,"Ciliate")
# diatom <- makeRatioPlot(dia,"Diatom")
# dinoflag <- makeRatioPlot(dino,"Dinoflagellate")
# haptophy <- makeRatioPlot(hapto,"Haptophyte")
# rhizaria <- makeRatioPlot(rhiz,"Rhizaria")

# Combine plots together as panels, create a common legend, and add panel labels (a-g)
allEuk+ (chlorophyte+ dinoflag+rhizaria)/(haptophy+ciliate+diatom)+plot_layout(guides = "collect",widths=c(2,4))+plot_annotation(tag_levels="a")
# ggsave("Figure5.pdf",width=26,height=11)


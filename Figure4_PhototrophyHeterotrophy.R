### Figure 4 - Phototrophy/Heterotrophy plots ###
### This script will allow you to summarize and plot the average +/- SE of various phototrophy/heterotrophy biomarker gene groups given a set of normalized transcript abundances ###
### Last Updated: July 24, 2023 ###

# Load libraries 
library(ggplot2)
library(tidyverse)
library(reshape2)
library(patchwork)

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
####

# Melt dataframe to make it long format
dfMelt <- melt(dfSum,id.vars="KO")

# Create new columns with eddy type and depth
cols <- colsplit(dfMelt$variable,"\\.",c("X","Eddy","Depth","Rep"))
dfMelt$Eddy <- cols$Eddy
dfMelt$Depth <- cols$Depth

# Gene groups - This is a list of KO terms that were defined as "phototrophy-associated" and "heterotrophy-associated" in this study
ko <- read.csv("KEGG.csv",header=TRUE)

# Join transcript abundance data and KO term list
dfMelt <- left_join(dfMelt,ko)

# Depth plot function
# Input: Dataframe containing transcript abundances and KO term list info, target gene group (column in KO list dataframe), and title of plot
# Output: Plot of mean Z-score transformed gene group abundance +/- SE across all eddy type/depth combinations
depthPlot <- function(df,target,title){
  df <- subset(df,Biomarker.Gene.Group==paste(target))
  df$KEGG.KO <- NULL
  df <- df[c(1:24)]
  dfMelt <- melt(df)
  dfSum <- dfMelt %>% group_by(variable) %>% summarize(sumz=sum(value)) %>% as.data.frame()
  m <- mean(dfSum$sumz)
  s <-sd(dfSum$sumz)
  dfZ <- dfSum %>% group_by(variable) %>% summarize(z=(sumz-m)/s)
  colz <- colsplit(dfZ$variable,"\\.",c("x","eddy","depth","rep"))
  dfZ$eddy <- colz$eddy
  dfZ$depth <- colz$depth
  
  dfFin <- dfZ %>% group_by(eddy,depth) %>% summarize(meanZ=mean(z),sdZ=(sd(z)/sqrt(3)))
  dfFin$depth2 <- c(122,150,250,25,112,150,250,25)
  
  dodge <- position_dodge(10)
  
  plot<- ggplot(dfFin, aes(x=depth2,y=meanZ,group=eddy,color=eddy))+geom_point(position=dodge,size=3.2)+scale_x_continuous(trans="reverse")+scale_y_continuous(position="right",limits=c(-3.5,3.5))+geom_errorbar(aes(ymin=meanZ-sdZ, ymax=meanZ+sdZ), width=0.4,position=dodge,linewidth=0.8)+coord_flip()+theme_classic()+scale_color_manual(name="Eddy",values=c("red","blue"))+xlab("Depth (m)")+ylab("Mean Z-score CPM +/- SE")+ggtitle(paste(title))+geom_hline(yintercept = 0,linetype="dotted")
  return(plot)}

# Implement function for the 12 gene groups focused on in this study
fat <- depthPlot(df,"Fatty acid breakdown","Fatty acid breakdown")
calvin <- depthPlot(df,"Calvin cycle","Calvin cycle")
n <- depthPlot(df,"Inorganic N uptake and assimilation","Inorganic N uptake + assimilation")
endo <- depthPlot(df,"Endocytosis","Endocytosis")
auto <- depthPlot(df,"Autophagy","Autophagy")
lyso <- depthPlot(df,"Lysosome","Lysosome")
cath <- depthPlot(df,"Cathepsin","Cathepsin Protease")
vtype <- depthPlot(df,"V-type ATPase","V-type ATPase")
photo <- depthPlot(df,"Photosynthesis","Photosynthesis")
phago <- depthPlot(df,"Phagosome","Phagosome")
rho <- depthPlot(df,"Rho GTPase","Rho GTPase")
ant <- depthPlot(df,"Antenna proteins","Antenna proteins")

# Combine all plots together with a common legend and add panel labels (a-l)
photo+calvin+n+ant+phago+vtype+fat+cath+rho+endo+lyso+auto+plot_layout(guides = "collect",nrow=3)+plot_annotation(tag_levels="a")
# ggsave("Figure4.pdf",width=15,height=8)

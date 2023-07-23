### Figure 2 - Taxa barplot ###
### This script will allow you to make barplots that show the relative abundances of transcripts associated with different taxonomic groups ###
### Last Updated: July 23, 2023 ###

# Load libraries 
library(ggplot2)
library(patchwork)
library(randomcoloR)
library(tidyverse)
library(reshape2)

# MORE TO COME - SHOW TAX ASSIGNMENTS AND AVERAGING #
# tax <- read.csv("./Figure2/TaxaBarplot_MAY2023.csv",header=TRUE,row.names=1)
# tax <- melt(tax,id.vars="fin")
# namez <- colsplit(tax$variable,"\\.",c("Num","Eddy","Depth","Rep"))
# tax$Eddy <- namez$Eddy
# tax$Depth <- namez$Depth
# tax$EddyDepth <- paste(tax$Eddy,tax$Depth,sep="_")
# tax$Eddy <- NULL
# tax$Depth <- NULL
# taxMean <- tax %>% group_by(EddyDepth,fin)%>%summarize(m=mean(value))%>%as.data.frame()

# Identify groups and group labels (eddy type, depth)
groups <- c("Cyclonic_250m", "Anticyclonic_250m", "Cyclonic_150m", "Anticyclonic_150m","Cyclonic_112m","Anticyclonic_112m","Cyclonic_25m","Anticyclonic_25m")
group_labels<-c("Cyclonic 250 m", "Anticyclonic 250 m", "Cyclonic 150 m", "Anticyclonic 150 m", "Cyclonic DCM (112 m)","Anticyclonic DCM (122 m)", "Cyclonic 25 m", "Anticyclonic 25 m")
taxMean$sample<-factor(taxMean$EddyDepth, levels=groups, labels = group_labels)

# Subset cyclonic and anticyclonic eddy samples
plotCy <- subset(taxMean, grepl("Cyclon",taxMean$sample))
plotAcy <- subset(taxMean, grepl("Anticyclon",taxMean$sample))

# Plot function
# Input: Dataframe, plot tite
# Output: Barplot
colrs <- randomcoloR::distinctColorPalette(length(unique(tax$fin)))
barplotFxn <- function(df,title){
  ggplot(df,aes(y=m,x=sample,fill=fin))+
    geom_bar(stat="identity", position="fill", color="#525252")+
    labs(title="", x="",y="Relative transcript abundance")+
    scale_x_discrete(limits=c(), expand = c(0, 0))+
    scale_fill_manual(name="Taxonomic Group", values=colrs)+
    coord_flip()+
    theme(legend.title=element_blank(),legend.position="right",legend.text.align=0, axis.text = element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(), axis.line = element_line())+theme_classic(base_size=14)+ggtitle(paste(title))}

# Apply plot function to cyclonic and anticyclonic eddy samples 
cy <- barplotFxn(plotCy,"Cyclonic")
acy <- barplotFxn(plotAcy,"Anticyclonic")

# Combine plots
cy+acy+plot_layout(guides = "collect",nrow=2)+plot_annotation(tag_levels="a")
#ggsave("Figure2.pdf",width=8,height=5)

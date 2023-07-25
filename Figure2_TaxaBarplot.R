### Figure 2 - Taxa barplot ###
### This script will allow you to make barplots that show the relative abundances of transcripts associated with different taxonomic groups ###
### Last Updated: July 24, 2023 ###

# Load libraries 
library(ggplot2)
library(patchwork)
library(randomcoloR)
library(tidyverse)
library(reshape2)

# Load in normalized dataframe (obtained via Normalize_Data.R)
df <- read.csv("normalized_metaT_data.csv",header=TRUE,row.names=1)

# Assign taxonomic groups of interest
df$tax <- ifelse(grepl("Chlorophyta",df$Taxonomy),"Chlorophyte",NA)
df$tax <- ifelse(grepl("Ciliophora",df$Taxonomy),"Ciliate",df$tax) 
df$tax <- ifelse(grepl("Cryptophyta",df$Taxonomy),"Cryptophyte",df$tax)
df$tax <- ifelse(grepl("Bacillariophyta",df$Taxonomy),"Diatom",df$tax)
df$tax <- ifelse(grepl("Dinophyceae",df$Taxonomy),"Dinoflagellate",df$tax)  
df$tax <- ifelse(grepl("Haptophyta",df$Taxonomy),"Haptophyte",df$tax)  
df$tax <- ifelse(grepl("Pelagophyceae",df$Taxonomy),"Pelagophyte",df$tax)
df$tax <- ifelse(grepl("Rhizaria",df$Taxonomy),"Rhizaria",df$tax)        
df$tax <- ifelse(grepl("Syndiniales",df$Taxonomy),"Syndiniales",df$tax)
df$tax <- ifelse(grepl("Opisthokont",df$Taxonomy),"Opisthokont",df$tax)   
df$tax <- ifelse(grepl("Alveolata",df$Taxonomy) & is.na(df$tax),"Other Alveolate",df$tax)
df$tax <- ifelse(grepl("Archaeplastida",df$Taxonomy) & is.na(df$tax),"Other Archaeplastida",df$tax)
df$tax <- ifelse(grepl("Stramenopiles",df$Taxonomy) & is.na(df$tax),"Other Stramenopiles",df$tax)  
df$tax <- ifelse(df$Taxonomy=="Eukaryota","Unknown Eukaryote",df$tax)
df$tax <- ifelse(is.na(df$tax),"Other Eukaryote",df$tax)

# Remove columns we do not need for taxonomy summary
df$contigID <- NULL
df$Taxonomy <- NULL
df$Cluster <- NULL
df$KO <- NULL

# Sum the number of transcripts per taxonomic group per sample
dfSum <- df %>% group_by(tax) %>% summarize_all(sum) %>% as.data.frame()

# Melt the dataframe to make it long
dfMelt <- melt(dfSum,id.vars="tax")

# Create new columns with eddy type and depth
cols <- colsplit(dfMelt$variable,"\\.",c("X","Eddy","Depth","Rep"))
dfMelt$Eddy <- cols$Eddy
dfMelt$Depth <- cols$Depth

# Calculate average abundances across reps for each eddy/depth
dfAvg <- dfMelt %>% group_by(tax,Eddy,Depth) %>% summarize(m=mean(value)) %>% as.data.frame()

# Create sample column
dfAvg$Sample <- paste(dfAvg$Eddy,dfAvg$Depth,sep="_")

# Remove columns that are not needed
dfAvg$Eddy <- NULL
dfAvg$Depth <- NULL

# Identify groups and group labels (eddy type, depth)
groups <- c("Cyclonic_250m", "Anticyclonic_250m", "Cyclonic_150m", "Anticyclonic_150m","Cyclonic_112m","Anticyclonic_122m","Cyclonic_25m","Anticyclonic_25m")
group_labels<-c("Cyclonic 250 m", "Anticyclonic 250 m", "Cyclonic 150 m", "Anticyclonic 150 m", "Cyclonic DCM (112 m)","Anticyclonic DCM (122 m)", "Cyclonic 25 m", "Anticyclonic 25 m")
dfAvg$Sample<-factor(dfAvg$Sample, levels=groups, labels = group_labels)

# Subset cyclonic and anticyclonic eddy samples
plotCy <- subset(dfAvg, grepl("Cyclon",dfAvg$Sample))
plotAcy <- subset(dfAvg, grepl("Anticyclon",dfAvg$Sample))

# Plot function
# Input: Dataframe, plot tite
# Output: Barplot
colrs <- randomcoloR::distinctColorPalette(length(unique(dfAvg$tax)))
barplotFxn <- function(df,title){
  ggplot(df,aes(y=m,x=Sample,fill=tax))+
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
# ggsave("Figure2.pdf",width=8,height=5)

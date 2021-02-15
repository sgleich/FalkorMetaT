###### Data Processing ######
# By: Samantha Gleich
# Parts of code modified from Sarah K. Hu (https://github.com/shu251/SPOT_metatranscriptome)
# Contributors: Sarah K. Hu

# Libraries needed
library(tidverse)
library(edgeR)
library(reshape2)

####### 1.) Compile all blastx outputs (taxonomy), salmon outputs (transcript abundance counts), ghostKO outputs (KEGG functional annotation), and eggNOG outputs (eggNOG functional annotation) ######

# These loops will take individual outputs from blastx, salmon, ghostKO, and eggNOG and will compile them into final files (blastx_final.csv, salmon_final.csv, ghostko_final.csv, eggnog_final.csv)

# Make empty compiled output df - blastx
blastx_final <- NULL

# List groups of interest (filenames)
groups <- c("anticyclonic_25","anticyclonic_112","anticyclonic_150","anticyclonic_250","cyclonic_25","cyclonic_112","cyclonic_150","cyclonic_250")

# For loop to combine blastx files with sample # info
for (group in groups){
  blastx_in<-read.delim(paste(group,"tax.txt",sep="_"), header = TRUE,row.names=NULL)
  colnames(blastx_in) <- c("contigID", "Taxonomy")
  blastx_in$Sample <- group
  blastx_final <- rbind(blastx_final, blastx_in)}

# write.csv(blastx_final,"blastx_final.csv")

# List samples (filenames)
samples <- c("10-Cyclonic-112m-A","11-Cyclonic-112m-B","12-Cyclonic-112m-C","13-Cyclonic-25m-A","14-Cyclonic-25m-B","15-Cyclonic-25m-C","19-Anticyclonic-250m-A","20-Anticyclonic-250m-B","21-Anticyclonic-250m-C","22-Anticyclonic-150m-A","23-Anticyclonic-150m-B","24-Anticyclonic-150m-C","25-Anticyclonic-112m-A","26-Anticyclonic-112m-B","27-Anticyclonic-112m-C","28-Anticyclonic-25m-A","29-Anticyclonic-25m-B","30-Anticyclonic-25m-C","4-Cyclonic-250m-A","5-Cyclonic-250m-B","6-Cyclonic-250m-C","7-Cyclonic-150m-A","8-Cyclonic-150m-B","9-Cyclonic-150m-C")

# Make empty compiled output df - salmon
salmon_compile <- NULL

# For loop to combine count files with sample # info
for (sample in samples){
  salmon_in<-read.delim(paste(sample,"_quant.sf",sep="_"), header = TRUE)
  salmon_in$SampleName <- sample
  salmon_compile <- rbind(salmon_compile, salmon_in)}

# write.csv(salmon_compile,"salmon_final.csv")


# Make empty compiled output df - ghostKO
ghostko_final <- NULL

# List the groups of interest (filenames)
groups <- c("anticyclonic_25","anticyclonic_112","anticyclonic_150","anticyclonic_250","cyclonic_25","cyclonic_112","cyclonic_150","cyclonic_250")

# For loop to combine ghostKO files with sample # info
for (group in groups){
  ko_in<-read.delim(paste(group,"ko_final.txt",sep="_"), header = TRUE,row.names=NULL)
  colnames(ko_in) <- c("contigID", "KO")
  ghostko_final <- rbind(ghostko_final, ko_in)}

# Save compiled output file
# write.csv(ghostko_final,"ghostko_final.csv")

# Make empty compiled output df - eggNOG
eggnog_final <- NULL

# List the groups of interest (filenames)
groups <- c("anticyclonic_25","anticyclonic_112","anticyclonic_150","anticyclonic_250","cyclonic_25","cyclonic_112","cyclonic_150","cyclonic_250")

# For loop to combine ghostKO files with sample # info
for (group in groups){
  nog_in<-read.delim(paste(group,"eggnog_annot.txt",sep="_"), header = TRUE,row.names=NULL)
  nog_in <- data.frame(nog_in$V1, nog_in$V19)
  best_og <- sub('.*\\,', '', nog_in$nog_in.V19)
  nog_in$OG <- best_og
  nog_in$nog_in.V19 <- NULL
  colnames(nog_in) <- c("contigID", "OG")
  eggnog_final <- rbind(eggnog_final, nog_in)}

# Save compiled output file
# write.csv(eggnog_final,"eggnog_final.csv")

###### 2.) Compile salmon, blastx, ghostKO, and eggNOG outputs into long and wide format dataframe ######

# These dataframes will be used in downstream analyses (e.g. edgeR) 

# blastx_final <- read.csv("blastx_final.csv", header=TRUE,row.names=1)
# ghostko_final <- read.csv("ghostko_final.csv", header=TRUE,row.names=1)
# eggnog_final <- read.csv("eggnog_final.csv", header=TRUE,row.names=1)
# salmon_final <- read.csv("salmon_final.csv", header=TRUE,row.names=1)

# Join taxonomy, KEGG IDs, and eggNOG OGs
tax.ko <- full_join(blastx_final,ghostko_final)
tax.ko <- left_join (tax.ko,eggnog_final) 

# Here, we will use eggNOG OGs only if there is no KEGG ID assigned to a specific contig
tax.ko$KO <- as.character(tax.ko$KO)
tax.ko$OG <- as.character(tax.ko$OG)
tax.ko$KO <- ifelse(is.na(tax.ko$KO), tax.ko$OG, tax.ko$KO)
tax.ko$OG <- NULL

# Add salmon counts to dataframe
# Make long dataframe first, then make wide format dataframe
salmon_final <- data.frame(salmon_final$Name,salmon_final$NumReads,salmon_final$SampleName)
colnames(salmon_final) <- c("contigID","NumReads","SampleName")
tax.ko.count <- left_join(tax.ko, salmon_final)
tax.ko.count <- subset(tax.ko.count, NumReads !=0)

# write.csv(tax.ko.count, "long_format_data.csv")

# Wide format dataframe
wide.data <- dcast(tax.ko.count, contigID+Taxonomy+KO~SampleName,fun.aggregate = sum,value.var= "NumReads")

# write.csv(wide.data,"wide_format_data.csv")

###### 3.) Add cluster info to wide data frame (python script merged_clust_w_df.py) ######

###### 4.) Normalize the data in edgeR to get CPM ######
# First, we will normalize the transcript abundance data for the whole dataset. Then, we will normalize the transcript abundance data for individual taxonomic groups we care about
# We will use the wide_format dataframe when doing edgeR normalization

# df <- read.csv("wide_format_data.csv", header=TRUE, row.names=1)

# Label columns 
colnames(df) <- c("contigID","Taxonomy","KO","Cluster","10-Cyclonic-112m-A","11-Cyclonic-112m-B","12-Cyclonic-112m-C","13-Cyclonic-25m-A","14-Cyclonic-25m-B","15-Cyclonic-25m-C","19-Anticyclonic-250m-A","20-Anticyclonic-250m-B","21-Anticyclonic-250m-C","22-Anticyclonic-150m-A","23-Anticyclonic-150m-B","24-Anticyclonic-150m-C","25-Anticyclonic-112m-A","26-Anticyclonic-112m-B","27-Anticyclonic-112m-C","28-Anticyclonic-25m-A","29-Anticyclonic-25m-B","30-Anticyclonic-25m-C","4-Cyclonic-250m-A","5-Cyclonic-250m-B","6-Cyclonic-250m-C","7-Cyclonic-150m-A","8-Cyclonic-150m-B","9-Cyclonic-150m-C")

# "Genes" = contigID, Taxonomy, KO, and Cluster
# "Groups" = groups of interest; here we have 8 eddy and depth combinations
dge.metaT <- DGEList(counts=df[5:28],genes=df[1:4],group=c(rep("Cyclonic-DCM",3),rep("Cyclonic-25m",3),rep("Anticyclonic-250m",3),rep("Anticyclonic-150m",3),rep("Anticyclonic-DCM",3),rep("Anticyclonic-25m",3),rep("Cyclonic-250m",3),rep("Cyclonic-150m",3)))

# dge.metaT$samples

# Use TMM normalization method
metaT.norm <- calcNormFactors(dge.metaT,method = "TMM")
cpm_data <- cpm(metaT.norm, normalized.lib.sizes=TRUE, log=FALSE) 
cpm_data<-as.data.frame(cpm_data)  
data_CPM<-data.frame(metaT.norm$genes,cpm_data)

# Calculate mean transcript CPM across eddy/depth replicates
data_CPM$MeanCPM_CyclonicDCM<-apply(data_CPM[5:7],1,FUN=mean)
data_CPM$MeanCPM_Cyclonic25<-apply(data_CPM[8:10],1,FUN=mean)
data_CPM$MeanCPM_Anticyclonic250<-apply(data_CPM[11:13],1,FUN=mean)
data_CPM$MeanCPM_Anticyclonic150<-apply(data_CPM[14:16],1,FUN=mean)
data_CPM$MeanCPM_AnticyclonicDCM<-apply(data_CPM[17:19],1,FUN=mean)
data_CPM$MeanCPM_Anticyclonic25<-apply(data_CPM[20:22],1,FUN=mean)
data_CPM$MeanCPM_Cyclonic250<-apply(data_CPM[23:25],1,FUN=mean)
data_CPM$MeanCPM_Cyclonic150<-apply(data_CPM[26:28],1,FUN=mean)

# write.csv(data_CPM, "data_CPM.csv")

###### 5.) Normalize the data for each taxonomic group of interest ######
# You can put this section of code in a script, change the taxonomic group of interest, change the file output name, and then rerun this script for all of the different taxonomic groups you're interested in

df <- read.delim("wide_format_data.csv", header=FALSE,sep="\t")

colnames(df) <- c("contigID","Taxonomy","KO","Cluster","10-Cyclonic-112m-A","11-Cyclonic-112m-B","12-Cyclonic-112m-C","13-Cyclonic-25m-A","14-Cyclonic-25m-B","15-Cyclonic-25m-C","19-Anticyclonic-250m-A","20-Anticyclonic-250m-B","21-Anticyclonic-250m-C","22-Anticyclonic-150m-A","23-Anticyclonic-150m-B","24-Anticyclonic-150m-C","25-Anticyclonic-112m-A","26-Anticyclonic-112m-B","27-Anticyclonic-112m-C","28-Anticyclonic-25m-A","29-Anticyclonic-25m-B","30-Anticyclonic-25m-C","4-Cyclonic-250m-A","5-Cyclonic-250m-B","6-Cyclonic-250m-C","7-Cyclonic-150m-A","8-Cyclonic-150m-B","9-Cyclonic-150m-C")

# Break down long "Taxonomy" column into individual columns
taxtmp <- colsplit(df$Taxonomy, ";", c("Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
df.tax <- cbind(df,taxtmp)

# Subset dataframe so that we only have the taxa of interest. Here we will subset for diatoms (Bacillariophyceae)
taxa.norm <- subset(df.tax, Class=="Bacillariophyceae")

# Normalize diatom-only dataframe (i.e. taxa.norm)
dge.metaT <- DGEList(counts=taxa.norm[5:28],genes=taxa.norm[1:4],group=c(rep("Cyclonic-DCM",3),rep("Cyclonic-25m",3),rep("Anticyclonic-250m",3),rep("Anticyclonic-150m",3),rep("Anticyclonic-DCM",3),rep("Anticyclonic-25m",3),rep("Cyclonic-250m",3),rep("Cyclonic-150m",3)))

metaT.norm <- calcNormFactors(dge.metaT,method = "TMM")
cpm_data <- cpm (metaT.norm, normalized.lib.sizes=TRUE, log=FALSE) 
cpm_data<-as.data.frame(cpm_data)  
data_CPM<-data.frame(metaT.norm$genes,cpm_data)

# Calculate mean CPM across replicates
data_CPM$MeanCPM_CyclonicDCM<-apply(data_CPM[5:7],1,FUN=mean)
data_CPM$MeanCPM_Cyclonic25<-apply(data_CPM[8:10],1,FUN=mean)
data_CPM$MeanCPM_Anticyclonic250<-apply(data_CPM[11:13],1,FUN=mean)
data_CPM$MeanCPM_Anticyclonic150<-apply(data_CPM[14:16],1,FUN=mean)
data_CPM$MeanCPM_AnticyclonicDCM<-apply(data_CPM[17:19],1,FUN=mean)
data_CPM$MeanCPM_Anticyclonic25<-apply(data_CPM[20:22],1,FUN=mean)
data_CPM$MeanCPM_Cyclonic250<-apply(data_CPM[23:25],1,FUN=mean)
data_CPM$MeanCPM_Cyclonic150<-apply(data_CPM[26:28],1,FUN=mean)

# write.csv(data_CPM, "diatom_data_CPM.csv")



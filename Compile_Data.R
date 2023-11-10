### This script will allow you to combine the salmon (count data), EUKulele (taxonomy data), eggNOG (function/gene data), and ortholog group (MMSeqs data) metaT outputs into one wide dataframe. This wide dataframe will be noramlized via (NormalizeData.R) ###
### By: Samantha Gleich ###
### Last Updated: November 10, 2023 ###

# Load libraries
library(tidyverse)
library(reshape2)

# Load in salmon data
salmon_fin <- read.csv("salmon_final.csv",header=TRUE,row.names=1)
salmon_fin2 <- data.frame(salmon_fin$Name, salmon_fin$NumReads,salmon_fin$SampleName)
colnames(salmon_fin2)<- c("Name", "NumReads","SampleName")

# Load in EUKulele taxonomy data
tax_fin <- read.csv("EUKulele_compiled.csv",header=TRUE,row.names=1)
tax_fin2 <- data.frame(tax_fin$Name, tax_fin$Taxonomy)
colnames(tax_fin2) <- c("Name","Taxonomy")

# Load in gene (i.e., functional annotation) information 
eggnog_fin <-read.csv("eggnog_final.csv",header=TRUE,row.names=1)

# Join eukulele, eggnog, and salmon outputs together. We can "join" the dataframes by the contig Name.
tax.ko <- full_join(tax_fin2,eggnog_fin)
tax.ko.count <- left_join(tax.ko, salmon_fin2)
all <- subset(tax.ko.count, NumReads !=0)

# Now we can save our data in long format (but we'll really want the data in wide format)
# write.csv(all, "long_format_data.csv")

# Let's wrangle the data so that it's in wide format
wide.data <- dcast(all, Name+Taxonomy+KEGG~SampleName,fun.aggregate = sum,value.var= "NumReads")

# Save wide format data
write.csv(wide.data,"wide_format_data.csv")

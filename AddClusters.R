### This script will allow you to combine the salmon (count data), EUKulele (taxonomy data), and eggNOG (function/gene data) metaT outputs into one wide dataframe. ###
### By: Samantha Gleich ###
### Last Updated: November 9, 2022 ###

# Load libraries
library(tidyverse)

# First let's load in our cluster file and make it more interpretable
df <- read.delim("clusters.tsv",header=FALSE,row.names=NULL)
df$V1 <- as.numeric(df$V1)
df2 <- df[c(2,1)]

# Save new cluster file (column #1 is contigID and column #2 is cluster #)
write.csv(df2, "clustFINAL.csv")

# Now we can load in our wide format data and join the cluster information to our dataframe with KO, taxonomy, and count information.
df1 <- read.csv("wide_format_data.csv",header=TRUE,row.names=1)
colnames(df2)<- c("contigID","Cluster")
df3 <- left_join(df1,df2)

# Reorder columns so that the cluster information is on the left-hand side by the contigID, KO information, and taxonomy information.
dffin <- df3[c(1:3,28,4:27)]
write.csv(dffin, "clustered_output_FINAL.csv")

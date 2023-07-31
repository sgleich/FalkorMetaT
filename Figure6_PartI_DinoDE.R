### Figure 6; Part I- Differential expression analyses ###
### This script will allow you to perform pairwise differentialy expression analyses in EdgeR ###
### Last Updated: July 31, 2023 ###

# Load libraries
library(edgeR) 

# Load in data
# First 2 columns in this dataframe are: orthologous group ID and KO annotation. Then the remaining columns (#3-26) are counts. Columns for counts correspond to different eddy types - depths - replicates (2 eddies * 4 depths * 3 reps = 24 columns of counts). 
df <- read.csv("Dino_Wide_Form.csv",header=TRUE,row.names=1)

# Set up differentially expressed gene list
dgList <- DGEList(counts=df[3:26],genes=df[1:2],group=c(rep("Cyclonic-DCM",3),rep("Cyclonic-25m",3),rep("Anticyclonic-250m",3),rep("Anticyclonic-150m",3),rep("Anticyclonic-DCM",3),rep("Anticyclonic-25m",3),rep("Cyclonic-250m",3),rep("Cyclonic-150m",3)))
dgList <- calcNormFactors(dgList,method = "TMM") # Calculate library normalization factors to scale library sizes
design <- model.matrix(~0+group, data=dgList$samples) # Means model
colnames(design) <- levels(dgList$samples$group)

# Differential expression model 
dgList <- estimateGLMCommonDisp(dgList,design) 
dgList <- estimateGLMTrendedDisp(dgList,design)
dgList <- estimateGLMTagwiseDisp(dgList,design)
fit <- glmQLFit(dgList, design)

# Set up pairwise comparisons. For the sake of Fig 6, we will look at cyclonic (1) vs. anticyclonic (-1) DGE at each depth. 
qlf <- glmQLFTest(fit, contrast=c(0,0,0,-1,0,0,0,1)) # Zeros = depths not being analyzed = model coefficients will be 0 for these terms.
qlf$comparison # Look at the comparison (e.g., -1*Anticyclonic-150m 1*Cyclonic-150m) 

# Obtain table of differentially expressed genes
table <- data.frame(qlf$table)
gene<- data.frame(qlf$genes)
total <- cbind(table,gene)
total$padj <- p.adjust(total$PValue,"fdr") # Adjust p-values

# Save output
write.csv(total,"DE_150m_Falkor.csv")


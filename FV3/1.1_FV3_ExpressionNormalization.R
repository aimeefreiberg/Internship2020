# 1.1 Fütterungs Versuch 3 Data Extraction - Aimee Freiberg 9.2020

#load in expression and annotation data
setwd("/home/aimee/NAS/Mouse/RNA/FV3/Intermediate")
FV3expr <- read.table("summaryArrays.txt", sep="\t", header = TRUE, check.names = FALSE)

setwd("/home/aimee/FV3 Data")
mouse <- read.table("FV3_Mice.tsv", sep="\t", header = TRUE)

#load necesssary libraries
library("affyPLM")
library("preprocessCore")

#see if data is normalized 
boxpolot(FV3expr)

#normalization of expression data
expressions <- normalize.quantiles(as.matrix(FV3expr))
colnames(expressions) <- colnames(FV3expr)
rownames(expressions) <- rownames(FV3expr)
expressions[1:4,]

#write out normalized data

write.table(expressions, "FV3_NormExpr.txt", sep="\t", col.names = TRUE)
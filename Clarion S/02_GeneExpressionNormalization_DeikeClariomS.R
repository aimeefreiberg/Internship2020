### load in expresssion data ###

setwd("/home/aimee/Aimee Test Data")
expressions1 <- read.csv("Expressionlevels_DeikeClarionS.txt", sep = "\t", header = TRUE)
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character", row.names=1)


library("affyPLM")
library("preprocessCore")

#normalize data 
expressions <- normalize.quantiles(as.matrix(expressions1))
colnames(expressions) <- colnames(expressions1)
rownames(expressions) <- rownames(expressions1)
expressions[1:4,]

#extract different tissues
ids <- unlist(lapply(strsplit(colnames(expressions)[1:nrow(arraymapping)], "_"),"[",2))
colnames(expressions)[1:nrow(arraymapping)] <- arraymapping[ids, 1]
colnames(expressions)[(nrow(arraymapping)+1):ncol(expressions)] = paste0("B6_", 1:8, "_L")

write.table(expressions, "Expressionlevels_DeikeClarionS.norm.samples.txt", sep = "\t", quote = FALSE, row.names = TRUE)

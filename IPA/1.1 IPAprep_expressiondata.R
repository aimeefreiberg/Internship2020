#IPA expression data preparation

setwd("/home/aimee/Aimee Test Data")
expressions <- read.csv("Expressionlevels_DeikeClarionS.norm.samples.txt", sep = "\t", header = TRUE)
ids = colnames(expressions)

#extract liver expression data <- adjust this so any tissue can be used for that  (P - Pankreas , SM skeletal Mucle, GF - Gonadal fat, L - Liver)
tissue <- which(grepl("L", ids))

#test for significant differences in tissue expression data (S1/S2/B6)
tissueexpr = expressions[, tissue]
expressions = expressions[, -tissue]

tissueexpr[1:5,]

strain = substr(colnames(tissueexpr), 0,2)

results = t(apply(tissueexpr, 1, function(Y){
  pval = anova(lm(Y ~ strain))[[5]][1]
  return(c(coef(lm(Y ~ strain + 0)), pval))
}))

rownames(results) <- gsub("_at", "", rownames(results))
results[1:5,]

# load in table of fibrosis genes 
setwd("/home/aimee/Aimee Test Data")
fibgenes <- read.table("fibrosisgenes.txt", sep="\t", header=TRUE, row.names=1, fill=TRUE)

# subset for the fibrosis genes
fibresults = results[which(rownames(results) %in% fibgenes[, "ensembl_gene_id"]), ]
colnames(fibresults) = c("meanB6", "meanS1", "meanS2", "p.value")
fibresults[1:5, ]

#calculate fold change 
fc <- log2(fibresults[,"meanS1"]/fibresults[,"meanS2"])

#prep IPA Table fc + pvalue
ExpressionIPA <- cbind(fc, fibresults[,"p.value"])


#annotate file entrez gene ID
annotate <- function(significant){
  library(biomaRt)
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  res.biomart <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"), 
                       filters = c("ensembl_gene_id"), 
                       values = rownames(significant), mart = bio.mart)
  return(res.biomart)
}

ExpressionAnno = annotate(ExpressionIPA)
ExpressionAnno[1:5, ]

Annotation = ExpressionAnno[- which(duplicated(ExpressionAnno[,1])), ] 


#bind relvant data into table

# bind data into one table  

rownames(Annotation) = Annotation[, 1]
IPAexpression <- cbind(Annotation[, 2], ExpressionIPA[rownames(Annotation), ])
colnames(IPAexpression) = c("entrezID","foldchange","p.value")

#write out table
write.table(IPAexpression, file = "S1vsS2.results_IPAformat.txt", sep="\t", col.names = TRUE)

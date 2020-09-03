# subsetting for fibrosis genes and annotating
# load in data
setwd("/home/aimee/Aimee Test Data")
fibgenes <- read.table("fibrosisgenes.txt", sep="\t", header=TRUE, row.names=1, fill=TRUE)

setwd("/home/aimee/FV3 Data")
results <- read.table("FV3_ExpressionAnalizedAll.txt", sep="\t", header = TRUE )

setwd("/home/aimee/NAS/Mouse/RNA/FV3/Annotation")
annotation <- read.table("probeannotation.txt", sep="\t", header = TRUE)

#match probeID to ensemble ID
NewAnno = annotation[!duplicated(annotation[, "ProbeName"]), ]
rownames(NewAnno) = NewAnno[, "ProbeName"]
NewRes = cbind(NewAnno[rownames(results), c("ensembl_gene_id", "mgi_symbol", "ProbeChr", "ProbeStart", "ProbeStop")], results)
rownames(NewRes) = rownames(results)

#write out annotated table
setwd("/home/aimee/FV3 Data")
write.table(NewRes, "FV3_ExprAnnotated.txt", sep="\t", col.names = TRUE)
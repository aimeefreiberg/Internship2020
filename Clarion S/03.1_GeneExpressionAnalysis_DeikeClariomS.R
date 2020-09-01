# Gene analysis for liver expression data for B6,S1,S2 <- maybe should be changed so the loop runs for all tissues ? 

### load in normalized expresssion data ###

setwd("/home/aimee/Aimee Test Data")
expressions <- read.csv("Expressionlevels_DeikeClarionS.norm.samples.txt", sep = "\t", header = TRUE)
ids = colnames(expressions)

#extract liver expression data <- adjust this so any tissue can be used for that  (P - Pankreas , SM skeletal Mucle, GF - Gonadal fat)
tissue <- which(grepl("L", ids))

# make heatmat to see if correlation looks normal and no samples are switched (quality control switch)
corM <- cor(expressions)
png("heatmap.png", width = 2000, height = 2000)
heatmap(corM, scale = "none")
dev.off()

#test for significant differences in tissue expression data (S1/S2/B6)
tissueexpr = expressions[, tissue]
expressions = expressions[, -tissue]

tissueexpr[1:4,]

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

#annotate file with ensemble data
annotate <- function(significant){
  library(biomaRt)
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                       filters = c("ensembl_gene_id"), 
                       values = rownames(significant), mart = bio.mart)
  rownames(res.biomart) <- res.biomart[,1]
  annotated <- cbind(res.biomart[rownames(significant), ], significant)
  return(annotated)
}

fibresults = annotate(fibresults)
fibresults[1:5, ]

# specifically annotated data table for IPA 

<- 

# extract data table for the tissue -> change file name for specific tissue
write.table(fibresults, "ResultsLiver_Expression_DeikeClariomS.txt", sep="\t", row.names = TRUE)
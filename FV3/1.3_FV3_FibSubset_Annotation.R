# subsetting for fibrosis genes and annotating


###fibrosis genes 
# load in table of fibrosis genes 
setwd("/home/aimee/Aimee Test Data")
fibgenes <- read.table("fibrosisgenes.txt", sep="\t", header=TRUE, row.names=1, fill=TRUE)
setwd("/home/aimee/NAS/Mouse/RNA/FV3/Annotation")
annotation <- read.table("probeannotation.txt", sep="\t", header = TRUE)

# subset for the fibrosis genes
fibresults = results[which(rownames(results) %in% fibgenes[, "ensembl_gene_id"]), ]
colnames(fibresults) = c("meanB6", "meanS1", "meanS2", "p.value")
fibresults[1:5, ]

which(rownames(results) %in% annotation[, "ProbeName"])


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
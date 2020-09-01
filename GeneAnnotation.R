# get Ensemble gene ID
#8.07.20

#set work directive & load data
setwd("/home/aimee/Aimee Test Data")
genenames = read.table("Gene1.tsv", sep="\t", row.names= 1, header=TRUE, na.strings = c("","NA"))

#load libraries
library("biomaRt")

#set database and dataset to ensembl mouse gene data
ensembl <- useMart("ensembl", dataset ="mmusculus_gene_ensembl")

#load gene name data for selected genes

res.ensembl <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_symbol", "mgi_description", "entrezgene_id"), 
                        filters = c("mgi_symbol"), 
                        values = rownames(genenames), mart = ensembl)

# bind data into one table  
rownames(res.ensembl) = res.ensembl[, "mgi_symbol"]
GeneID <- cbind(genenames[,2:4], res.ensembl[, "ensembl_gene_id"][match(rownames(genenames), rownames(res.ensembl))])

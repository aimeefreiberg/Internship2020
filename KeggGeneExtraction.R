# get Genes from relevant Pathways for Fibrosis from KEGG -- 13.7.2020

# use bioconductor as source (only if needed to install for the first time)
# source("http://bioconductor.org/biocLite.R")

#paste the following code to install:
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")

#BiocManager::install("KEGGREST")

# load KEGGREST library
library(KEGGREST)

#find genes related to fibrosis in mice
#Fibgenes <- keggFind("mmu", "fibrosis")
#TGFb <- keggFind("mmu", "tgfb")

# pathways related to fibrosis in mice

fibpath <- c("mmu04932", "mmu05330", "mmu05410", "mmu04926", "mmu04614", "mmu04371")

# extract all the genes related to relevant pathways (make a foreloop)

fibrosisgenes = c()
for (path in fibpath) {
	genes <- keggLink("mmu",path)
	fibrosisgenes = c(fibrosisgenes, genes)
	}

mpathways = cbind(Pathway = names(fibrosisgenes), KEGG.ID = fibrosisgenes)
	
	
#convert KEGG ID to NGBI ID

NCBI.ID = c()
for (gene in mpathways[, 2]) {
	ID <- keggConv("ncbi-geneid", gene)
	NCBI.ID = c(NCBI.ID, ID)
	}
mpathways = cbind(mpathways, NCBI.ID)
mpathways[, 3] <- gsub("ncbi-geneid:", "", mpathways[,3])


#load libraries
library("biomaRt")

#set database and dataset to ensembl mouse gene data
ensembl <- useMart("ensembl", dataset ="mmusculus_gene_ensembl")


#load gene name data for selected genes

res.ensembl.new <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_symbol", "mgi_description", "entrezgene_id"), 
                        filters = c("entrezgene_id"), 
                        values = mpathways[, 3], mart = ensembl)


#look which filter is the ensembl ID & take the valid name
#grep("NCBI", listFilters(ensembl)[, 2])
#listFilters(ensembl)[c(26,81,82), ]

######## get Ensemble gene ID previous research######
#8.07.20

#set work directive & load data
setwd("/home/aimee/AimÃ©e Test Data/")
genenames = read.table("Gene1.tsv", sep="\t", row.names= 1, header=TRUE, na.strings = c("","NA"))

#load libraries
library("biomaRt")

#set database and dataset to ensembl mouse gene data
ensembl <- useMart("ensembl", dataset ="mmusculus_gene_ensembl")

#load gene name data for selected genes

res.ensembl <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_symbol", "mgi_description", "entrezgene_id"), 
                        filters = c("mgi_symbol"), 
                        values = rownames(genenames), mart = ensembl)

# bind all genes into one table
mydata <- rbind(res.ensembl, res.ensembl.new)

# write out gene table 
write.table(mydata, file = "fibrosisgenes.txt", sep = "\t")




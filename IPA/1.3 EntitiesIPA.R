# prepare Entities text file for IPA

#load in relevant filed
setwd("/home/aimee/IPA/Data")
Entities <- read.table("StringEnts.dat", sep="\t", header = TRUE )

#prep table for modification
colnames(Entities) = c("uid","id","ensembl_peptide_id","symbol")
Entities[,  "ensembl_peptide_id"] <- sapply(strsplit(Entities[, "ensembl_peptide_id"], split=".", fixed=TRUE), function(x) (x[2]))
rownames(Entities) = Entities[,  "ensembl_peptide_id" ]

# add entrez id and hgnc symbol matched to the ensemble protein ID

annotate <- function(significant){
  library(biomaRt)
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  res.biomart <- getBM(attributes = c("entrezgene_id","ensembl_peptide_id", "mgi_symbol"), 
                       filters = c("ensembl_peptide_id"), 
                       values = rownames(significant), mart = bio.mart)
  return(res.biomart)
}

AEntities<- annotate(Entities)
AEntities [1:5, ]

#insert id and symbol by matching ensembl gene ID

for (x in 1:nrow(Entities)){
			row <- which(AEntities[, "ensembl_peptide_id"] == Entities[x, "ensembl_peptide_id"])
			if(length(row)==1){
			Entities[x, "id"] <- AEntities[row, "entrezgene_id"]
			Entities[x, "symbol"] <- AEntities[row, "mgi_symbol"]
			}	
	}
	
Entities[1:5, ]

write.table(Entities, "EntitiesIPA.dat", sep="\t")
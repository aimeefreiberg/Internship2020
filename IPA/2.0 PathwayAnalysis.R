# Pathway analysis package QuarternaryProd alternative to IPA - Aimée Freiberg 8.2020


# download IPA package 
#if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#BiocManager::install("QuaternaryProd")

# load package 
library("QuaternaryProd")

# load in relevant data
setwd("/home/aimee/IPA/Data")
Expression <- read.table("S1vsS2.results_IPAformat_liver.txt", sep = "\t", header = TRUE)
Relations <- read.table("stringRels.dat", sep="\t", header = TRUE)
Entities <- read.table("EntitiesIPA.dat", sep="\t", header = TRUE)

#check that the entret is character or integral & p + fc is numeric
colnames(Expression) <- c("entrez", "pvalue", "fc")
colnames(Entities) <- c("uid", "id", "ensembleid", "symbol")

#remove duplicated rows from Relations Table
SubRelations <- Relations[!duplicated(Relations), ]

# remove NAs from Entitites

try <- apply(Entities, 1, is.na)
hasNA <- which(apply(try, 2, sum) > 0)
NAuid <- Entities[hasNA, "uid"]
Entities <- Entities[- hasNA, ]

#remove missing uid from Relations

uff <- c()
for (x in 1:nrow(SubRelations)) {
		if(SubRelations[x, "srcuid"] %in% NAuid){
		uff = c(uff, x)
		}
		if(SubRelations[x, "trguid"] %in% NAuid){
		uff = c(uff, x)
		}
}

SubRelations <- SubRelations[-unique(uff), ]

Entities = data.frame(Entities)
Entities[, "id"] = as.character(Entities[, "id"])
Entities[, "ensemblid"] = as.character(Entities[, "ensemblid"])
Entities[, "symbol"] = as.character(Entities[, "symbol"])

# run casual relation engine

QuartnernaryResults <- RunCRE_HSAStringDB(Expression,  
					 method = "Quaternary", 
					 fc.thresh = log2(1.3), 
					 pval.thresh = 0.05,
					 only.significant.pvalues = TRUE, 
					 significance.level = 0.05,
					 epsilon = 1e-16, 
					 progressBar = TRUE, 
					 relations = SubRelations, 
					 entities = Entities)

QuartnernaryResults[1:10, ]

#print histogram
png("histq.png")
hist(-log10(QuartnernaryResults[, "pvalue"]))
dev.off()

#subset for only significant pvalues and matching regulations

sig <- which(QuartnernaryResults[, "pvalue"] > 0 & QuartnernaryResults[, "pvalue"] < 0.05 & QuartnernaryResults[, "score"] == 1)
qres <- QuartnernaryResults[sig, ]

#write out significant regulators to match for pathways connected to fibrosis in innatedb.com
qres <- cbind(qres, Entities[qres[, "uid"], "id"])
write.table(qres, "regulators.txt", sep="\t")
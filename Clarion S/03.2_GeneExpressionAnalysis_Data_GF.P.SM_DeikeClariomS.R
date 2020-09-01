# analysie data for gonadal fat, skeletal muscle, and pancreas and select for the genes relevant for fibrosis

### load in normalized expresssion data ###

setwd("/home/aimee/Aimee Test Data")
expressions <- read.csv("Expressionlevels_DeikeClarionS.norm.samples.txt", sep = "\t", header = TRUE)
ids = colnames(expressions)

gonadalfat <- which(grepl("G", ids))
skeletalmuscle <- which(grepl("S", ids))
pankreas <- which(grepl("P", ids))

# make heatmat to see if corellation looks nomrla and no saples are switched (quality control switch)
corM <- cor(expressions)
png("heatmap.png", width = 2000, height = 2000)
heatmap(corM, scale = "none")
dev.off()

#get significant data when comparing both strains for the different tissues
getSignificant <- function(expressions, Tissue = "G", adjust = "BH", p.val = 0.05){
  S1P <- which(grepl("S1", colnames(expressions)) & grepl(Tissue, colnames(expressions)))
  S2P <- which(grepl("S2", colnames(expressions)) & grepl(Tissue, colnames(expressions)))

  res <- t(apply(expressions[, c(S1P, S2P)],1, function(x) {
    s1v <- x[1:length(S1P)]
    s2v <- x[(length(S1P)+1):(length(S1P) + length(S2P))]
    pval <- tryCatch({t.test(s1v, s2v)$p.value}, error = function(e) { return(1) })
    return(c(mean(s1v), mean(s2v), sd(s1v), sd(s2v), mean(s1v) / mean(s2v), log2(mean(s1v) / mean(s2v)), pval))
  }))
  colnames(res) <- c("mean(s1)", "mean(s2)", "sd(s1)", "sd(s2)", "FC", "logFC", "p.value")
  significant <- res[which(p.adjust(res[,"p.value"], adjust) < p.val),]
  #significant <- res
  rownames(significant) <- gsub("_at", "", rownames(significant))
  return(significant)
}

#extract pvalue for all the genes 
gonadalfat <- getSignificant(expressions, "G", p.val = 1.1)
skeletalmuscle <- getSignificant(expressions, "S", p.val = 1.1)
pankreas <- getSignificant(expressions, "P", p.val = 1.1)

#load in the fibrosis genes 
setwd("/home/aimee/Aimee Test Data")
fibgenes <- read.table("fibrosisgenes.txt", sep="\t", header=TRUE, row.names=1, fill=TRUE)

#subset fibrosis genes from the main table
fibGFat = gonadalfat[which(rownames(gonadalfat) %in% fibgenes[, "ensembl_gene_id"]), ]
fibSMuscle = skeletalmuscle[which(rownames(skeletalmuscle) %in% fibgenes[, "ensembl_gene_id"]), ]
fibpankreas = pankreas[which(rownames(pankreas)  %in% fibgenes[, "ensembl_gene_id"]), ]

#bind fibrosis gene data into one table
colnames(fibGFat)[colnames(fibGFat) %in% c("mean(s1)","mean(s2)","logFC","p.value")] = c("GFmean(s1)","GFmean(s2)","GFlogFC","GFp.value")
colnames(fibSMuscle)[colnames(fibSMuscle) %in% c("mean(s1)","mean(s2)","logFC","p.value")] = c("SMmean(s1)","SMmean(s2)","SMlogFC","SMp.value")
colnames(fibpankreas)[colnames(fibpankreas) %in% c("mean(s1)","mean(s2)","logFC","p.value")] = c("Pmean(s1)","Pmean(s2)","PlogFC","Pp.value")



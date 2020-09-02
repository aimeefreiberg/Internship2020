# Expression analysis Aimée Freiberg 2020

#load in normalized expression data 
setwd("/home/aimee/FV3 Data")
expr <- read.table("FV3_NormExpr.txt", sep="\t", header = TRUE)
mouse <- read.table("FV3_Mice.tsv", sep="\t", header = TRUE)

#extract for the different tissues (fat, liver, brain)
colnames(expr) <- sub("X", "", colnames(expr))
tissuemouse <- mouse[which(mouse[, "tissue"] == "fat"),]
tissueexpr <- expr[, gsub("-", "_", tissuemouse[, "core_name"])]

# test for significant difference in strains (B6 & BFMI860)
rownames(tissuemouse) = gsub("-", "_", tissuemouse[, "core_name"])
strain = tissuemouse[colnames(tissueexpr), "Direction.of.cross"]

results = t(apply(tissueexpr, 1, function(Y){
  pval = anova(lm(Y ~ strain))[[5]][1]
  return(c(coef(lm(Y ~ strain + 0)), pval))
}))

colnames(results) = c("B6", "BFRMI860", "pvalue")

write.table(results, "FV3_ExpressionAnalizedAll.txt", sep="\t", col.names = TRUE)
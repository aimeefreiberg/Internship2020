## work in progress 

####extract for the different tissues####

#fat 
fmouse <- mouse[which(mouse[, "tissue"] == "fat"),]
fatexpr <- expressions[, gsub("-", "_", fmouse[, "core_name"])]›
write.table(fatexpr, "FV3_fat.txt" sep)

#liver
lmouse <- mouse[which(mouse[, "tissue"] == "liver"),]

#brain
bmouse <- mouse[which(mouse[, "tissue"] == "brain"),]


write.table(expressions, "Expressionlevels_DeikeClarionS.norm.samples.txt", sep = "\t", quote = FALSE, row.names = TRUE)

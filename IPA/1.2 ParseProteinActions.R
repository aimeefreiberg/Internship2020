# Parse protein actions to StringRels.dat & StringEnts.dat
setwd("C:/Users/Arends/Downloads/10090.protein.actions.v11.0.txt")

mdata <- read.csv("10090.protein.actions.v11.0.txt",sep = "\t")

uniqueproteins <- unique(c(mdata[,"item_id_a"],mdata[,"item_id_b"]))

stringEnts <- cbind(uid = 1:length(uniqueproteins), id = NA, ensembleid = uniqueproteins, symbol = NA)

write.table(stringEnts, file = "StringEnts.dat", sep = "\t", quote = FALSE, row.names=FALSE)

msubset <- mdata[, c("item_id_a", "item_id_b", "mode")]
msubset <- msubset[-which(!msubset[,"mode"] %in% c("activation", "inhibition", "expression")),]

stringRels <- cbind(srcuid = unlist(lapply(msubset[,"item_id_a"], function(x){which(stringEnts[, "ensembleid"] == x)})), 
                    trguid = unlist(lapply(msubset[,"item_id_b"], function(x){which(stringEnts[, "ensembleid"] == x)})),
                    mode = unlist(lapply(msubset[, "mode"], function(x){if(x == "activation") return(1);  if(x == "inhibition") return(-1);  if(x == "expression") return(0);})))
write.table(stringRels, file = "stringRels.dat", sep = "\t", quote = FALSE,row.names=FALSE)

# extract the raw expression data from CEL files (Deike S1 S2 data + B6 Liver samples


library(affy)

setwd("/home/aimee/NAS/Mouse/RNA/ClariomS_S1vsS2_Deike/DN-2019_8745-Data/Rohdaten")

#read in the cel files for BFMI and B6N mie with annotation 
dat <- ReadAffy(cdfname ="clariomsmousemmensgcdf") # Use the clariomsmouse CDF from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/enst.asp

setwd("/home/aimee/B6_Liver")
bdat <- ReadAffy(cdfname ="clariomsmousemmensgcdf")

#call expression level based on annotation file
eset <- mas5(dat)
beset <- mas5(bdat)

#log 2 transform and extract the expression data
expr <- log2(assayData(eset)$exprs)
bexpr <- log2(assayData(beset)$exprs)

#bind them together in one table
alleset <- cbind(expr,bexpr)

alleset[1:4, ]

#write out expression data / read expression data

setwd("/home/aimee/Aimee Test Data")
write.table(alleset, "Expressionlevels_DeikeClarionS.txt", sep = "\t", quote = FALSE, row.names = TRUE)


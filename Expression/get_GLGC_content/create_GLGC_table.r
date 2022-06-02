library(EDASeq)
library(biomaRt)

#get ensembl gene IDs for hg38
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
biomart_getID <- getBM(attributes = c("ensembl_gene_id"), mart = ensembl)

#get gene length and GC content for all IDs
getdata <- getGeneLengthAndGCContent(biomart_getID$ensembl_gene_id , org="hsa", mode = c("biomart"))
save(getdata, file = "getGLGC_download.rda")

#Save output as data frame with correct header names
geneInfoHT <- data.frame(geneLength = getdata[,1] , 
                          gcContent = getdata[,2])
#Save final table
save(geneInfoHT, file = "geneInfoHT.rda")
write.table(geneInfoHT, file="GLGC_table.txt")
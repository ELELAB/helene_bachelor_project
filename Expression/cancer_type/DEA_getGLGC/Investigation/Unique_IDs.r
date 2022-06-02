
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

gdc <- get(load(file = "../../../../data/BRCA.exp.rda"))

old_geneInfoHT <- read.table(file = "../../../get_GLGC_content/oldGC_geneinfoHT.txt")
new_geneInfoHT <- get(load(file = "../../../get_GLGC_content/geneInfoHT.rda"))

#get ensembl id for the two tables
ID_old_geneInfoHT <- rownames(old_geneInfoHT)
ID_new_geneInfoHT <- rownames(new_geneInfoHT)

#320 genes are present in geneInfoHT but not in biomartGC
not_in_new <- setdiff(ID_old_geneInfoHT,ID_new_geneInfoHT)

#44839 genes are present in biomartGC and not in geneInfoHT
not_in_old <- setdiff(ID_new_geneInfoHT, ID_old_geneInfoHT)

#read files with unique IDs from the Venn diagram
dif_upreg_oldGCn <- readRDS("dif_upreg_oldGCn.RDS")
dif_downreg_oldGCn <- readRDS("dif_downreg_oldGCn.RDS")

dif_upreg <- readRDS("dif_upreg.RDS")
dif_downreg <- readRDS("dif_downreg.RDS")

#convert hugo to ensembl for all unique DEGs
dif_upreg_oldGCn_ensembl <- rowData(gdc)[match(dif_upreg_oldGCn, rowData(gdc)[,"external_gene_name"]),"ensembl_gene_id"]
dif_downreg_oldGCn_ensembl <- rowData(gdc)[match(dif_downreg_oldGCn, rowData(gdc)[,"external_gene_name"]),"ensembl_gene_id"]

dif_upreg_ensembl <- rowData(gdc)[match(dif_upreg, rowData(gdc)[,"external_gene_name"]),"ensembl_gene_id"]
dif_downreg_ensembl <- rowData(gdc)[match(dif_downreg, rowData(gdc)[,"external_gene_name"]),"ensembl_gene_id"]


#investigate if unique upregulated genes are present in both tables
unique_upreg_in_new <- na.omit(match(dif_upreg_ensembl, ID_old_geneInfoHT))
unique_upreg_in_old <- na.omit(match(dif_upreg_oldGCn_ensembl, ID_new_geneInfoHT))

#conclusion: 400 genes are unique upregulated in oldGCn and all are present in biomartGC
# 718 of 2434 unique new upregulated genes are present in geneInfoHT

unique_downreg_in_new <- na.omit(match(dif_downreg_ensembl, ID_old_geneInfoHT))
unique_downreg_in_old <- na.omit(match(dif_downreg_oldGCn_ensembl, ID_new_geneInfoHT))

#conclusion: all 105 genes unique downregulated in oldGCn are present in biomartGC
#1015 of 2267 unique new downregulated genes are present in geneInfoHT

#Collecting data of the unique IDs normalized by the getGLGC table
collection_getGLGC <- data.frame(DE = rep(c("Downreg","Upreg"),2),
                         Value = c(length(dif_downreg_ensembl),length(dif_upreg_ensembl),
                                   length(unique_downreg_in_new),length(unique_upreg_in_new)),
                         Category = c(rep("Total uniqe",2),rep("In geneInfoHT Table",2)))
#barplot                          
ggplot(collection_getGLGC, aes(x = DE, y =  Value, fill = Category)) + 
  geom_bar(position = "dodge", stat = "identity") +
  ggtitle("DEGs unique in Venndiagram but some present in geneInfoHT table")+
  theme(legend.text = element_text(size = 14),axis.text = element_text(size = 12), axis.title=element_text(size=14),legend.title=element_text(size=16),
        plot.title = element_text(size=16))
ggsave("Barplot_unique_IDs.png")

#Collecting data of the unique IDs normalized by the geneInfoHT table
collection_geneInfoHT <- data.frame(DE = rep(c("Downreg_oldGCn","Upreg_oldGCn"),2),
                                 Value = c(length(dif_downreg_oldGCn_ensembl),length(dif_upreg_oldGCn_ensembl),
                                           length(unique_downreg_in_old),length(unique_upreg_in_old)),
                                 Category = c(rep("Total uniqe",2),rep("In getGLGC Table",2)))

ggplot(collection_geneInfoHT, aes(x = DE, y = Value, fill = Category)) + 
  geom_bar(position = "dodge", stat = "identity") +
  ggtitle("oldGCn DEGs unique in Venndiagram but all present in getGLGC table")+
  theme(legend.text = element_text(size = 14),axis.text = element_text(size = 12), axis.title=element_text(size=14),legend.title=element_text(size=16),
        plot.title = element_text(size=16))
ggsave("Barplot_unique_IDs_oldGCn.png")


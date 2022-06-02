#This script compares differential expressed genes (DEGs), which only difference
#is their normalization.
#A Venn-diagram, correlation plots are created as well as a enrichment analysis

#load libraries
library(grid)
library(gridBase)
library(futile.logger)
library(VennDiagram)
library(enrichR)

#load files
DEG_ensembl <- get(load(file = "DEG_ensembl.rda" ))
DEG_ensembl_oldGCn <- get(load(file = "DEG_ensembl_oldGCn.rda" ))
DEG_HUGO <- get(load(file = "DEG_HUGO.rda" ))
DEG_HUGO_oldGCn <- get(load(file = "DEG_HUGO_oldGCn.rda" ))


#Sort data into up- and downregulated DEGs ENSEMBL

ensembl_pos <- DEG_ensembl[DEG_ensembl$logFC >0,]
ensembl_neg <- DEG_ensembl[DEG_ensembl$logFC<0,]

ensembl_oldGCn_pos <- DEG_ensembl_oldGCn[DEG_ensembl_oldGCn$logFC>0,]
ensembl_oldGCn_neg <- DEG_ensembl_oldGCn[DEG_ensembl_oldGCn$logFC<0,]

#Sort data into up- and downregulated DEGs HUGO

hugo_pos <- DEG_HUGO[DEG_HUGO$logFC >0,]
hugo_neg <- DEG_HUGO[DEG_HUGO$logFC<0,]

hugo_oldGCn_pos <- DEG_HUGO_oldGCn[DEG_HUGO_oldGCn$logFC>0,]
hugo_oldGCn_neg <- DEG_HUGO_oldGCn[DEG_HUGO_oldGCn$logFC<0,]

#Find hugo DEGs, which doesn't have a blank ID and save as txt files
hugo_pos_wID <- hugo_pos[which(hugo_pos$ID != ""),]
saveRDS(hugo_pos_wID , file = "Upreg_hugo.RDS")

hugo_neg_wID <- hugo_neg[which(hugo_neg$ID != ""),]
saveRDS(hugo_neg_wID, file = "Downreg_hugo.RDS" )

hugo_oldGCn_pos_wID <- hugo_oldGCn_pos[which(hugo_oldGCn_pos$ID != ""),]
saveRDS(hugo_oldGCn_pos_wID, file = "Upreg_hugo_oldGCn.RDS" )

hugo_oldGCn_neg_wID <-hugo_oldGCn_neg[which(hugo_oldGCn_neg$ID != ""),]
saveRDS(hugo_oldGCn_neg_wID, file = "Downreg_hugo_oldGCn.RDS" )


#Creating Venn diagrams
venn_ensembl <- venn.diagram(
  x = list(rownames(ensembl_pos), rownames(ensembl_neg),rownames(ensembl_oldGCn_pos),rownames(ensembl_oldGCn_neg)),
  category.names = c("Upreg" , "Downreg" , "Upreg oldGCn","Downreg oldGCn"),
  filename = NULL
  
)

venn_hugo <- venn.diagram(
  x = list(hugo_pos_wID$ID, hugo_neg_wID$ID, hugo_oldGCn_pos_wID$ID,hugo_oldGCn_neg_wID$ID ),
  category.names = c("Upreg" , "Downreg" , "Upreg oldGCn","Downreg oldGCn"),
  filename = NULL
  
)
#Save Venn diagrams
png('venn_ensembl.png', width = 12,height =8,res = 300, units = "cm")
grid.draw(venn_ensembl)
dev.off()

png('venn_hugo.png', width = 12,height =8,res = 300, units = "cm")
grid.draw(venn_hugo)
dev.off()


##########################################################################################################
#top 20 genes for up- and down regulation in ENSEMBL

ensembl_pos20 <-ensembl_pos[order(ensembl_pos$logFC,decreasing = TRUE), ] [1:20,]
ensembl_neg20 <- ensembl_neg[order(ensembl_neg$logFC), ] [1:20,]

ensembl_oldGCn_pos20 <-ensembl_oldGCn_pos[order(ensembl_oldGCn_pos$logFC,decreasing = TRUE), ] [1:20,]
ensembl_oldGCn_neg20 <- ensembl_oldGCn_neg[order(ensembl_oldGCn_neg$logFC), ] [1:20,]

#top 20 genes for up- and down regulation in HUGO

hugo_pos20 <-hugo_pos[order(hugo_pos$logFC,decreasing = TRUE), ] [1:20,]
hugo_neg20 <- hugo_neg[order(hugo_neg$logFC), ] [1:20,]

hugo_oldGCn_pos20 <-hugo_oldGCn_pos[order(hugo_oldGCn_pos$logFC,decreasing = TRUE), ] [1:20,]
hugo_oldGCn_neg20 <- hugo_oldGCn_neg[order(hugo_oldGCn_neg$logFC), ] [1:20,]

##########################################################################################################
#testing if its the same most up and downregulated genes in new and old analysis ENSEMBL

ensembl_compare_pos <- setequal(rownames(ensembl_pos20),rownames(ensembl_oldGCn_pos20))
ensembl_compare_neg <- setequal(rownames(ensembl_neg20),rownames(ensembl_oldGCn_neg20))

#testing if its the same most up and downregulated genes in new and old analysis HUGO

hugo_compare_pos <- setequal(rownames(hugo_pos20),rownames(hugo_oldGCn_pos20))
hugo_compare_neg <- setequal(rownames(hugo_neg20),rownames(hugo_oldGCn_neg20))

#Produce output of comparison results
if (ensembl_compare_pos == FALSE){
  ensembl_pos_string <- "The ensembl gene ID's upregulated are not identical between the two analysis"
}  else {
    ensembl_pos_string <- "The ensembl gene ID's upregulated are identical between the two analysis"
}

if (ensembl_compare_neg == FALSE){
  ensembl_neg_string <- "The ensembl gene ID's downregulated are not identical between the two analysis"
}  else {
  ensembl_neg_string <- "The ensembl gene ID's downregulated are identical between the two analysis"
}

if (hugo_compare_pos == FALSE){
  hugo_pos_string <- "The HUGO gene ID's upregulated are not identical between the two analysis"
}  else {
  hugo_pos_string <- "The HUGO gene ID's upregulated are identical between the two analysis"
}

if (hugo_compare_neg == FALSE){
  hugo_neg_string <- "The HUGO gene ID's downregulated are not identical between the two analysis"
}  else {
  hugo_neg_string <- "The HUGO gene ID's downregulated are identical between the two analysis"
} 

#save results in txt file

cat(ensembl_pos_string,ensembl_neg_string,hugo_pos_string,hugo_neg_string,sep="\n", file = "ID_comparison.txt")

##########################################################################################################

#Create correlation plot of ENSEMBL DEGs
png('correlationplot_ENSEMBL.png', width = 16,height =8,res = 300, units = "cm")
par(mfrow=c(1,2))
plot(ensembl_pos20$logFC, ensembl_oldGCn_pos20$logFC ,xlab = "New top 20 logFC",ylab = "Old top 20 logFC",main = "Upregulated ENSEMBL")
plot(ensembl_neg20$logFC, ensembl_oldGCn_neg20$logFC ,xlab = "New top 20 logFC",ylab = "Old top 20 logFC",main = "Downregulated ENSEMBL")
dev.off()

#Create correlation plot of ENSEMBL DEGs
png('correlationplot_HUGO.png', width = 16,height =8,res = 300, units = "cm")
par(mfrow=c(1,2))
plot(hugo_pos20$logFC, hugo_oldGCn_pos20$logFC ,xlab = "New top 20 logFC",ylab = "Old top 20 logFC",main = "Upregulated HUGO")
plot(hugo_neg20$logFC, hugo_oldGCn_neg20$logFC ,xlab = "New top 20 logFC",ylab = "Old top 20 logFC",main = "Downregulated HUGO")
dev.off()

##########################################################################################################
#Enrichment analysis
dbs <- c("GO_Biological_Process_2021","KEGG_2021_Human","Reactome_2016")

enriched_hugo_neg <- enrichr(hugo_neg_wID$ID , dbs)
enriched_hugo_pos <- enrichr(hugo_pos_wID$ID , dbs)
enriched_hugo_oldGCn_neg <- enrichr(hugo_oldGCn_neg_wID$ID , dbs)
enriched_hugo_oldGCn_pos <- enrichr(hugo_oldGCn_pos_wID$ID , dbs)

#Save enrichment plots

png('Enrichment_hugo_down.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_hugo_neg[["GO_Biological_Process_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_hugo_up.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_hugo_pos[["GO_Biological_Process_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_hugo_oldGCn_down.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_hugo_oldGCn_neg[["GO_Biological_Process_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_hugo_oldGCn_up.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_hugo_oldGCn_pos[["GO_Biological_Process_2021"]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()






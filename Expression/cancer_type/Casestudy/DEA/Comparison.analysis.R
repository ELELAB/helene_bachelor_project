#This script compares differential expressed genes (DEGs) obtained by
#two different methods: Limma-voom and edgeR GLM.
#A Venn-diagram, correlation plots are created as well as a enrichment analysis

#load libraries
library(grid)
library(gridBase)
library(futile.logger)
library(VennDiagram)
library(enrichR)

#load files and create ID column for edgeR DEGs
DEG_Limma_wID <- get(load(file = "DEG_Limma.rda" ))

DEG_edgeR_wID <- get(load(file = "DEG_edgeR.rda" ))
DEG_edgeR_wID$ID <- rownames(DEG_edgeR_wID) 

#Sort data into up- and downregulated DEGs and save files

Limma_pos <- DEG_Limma_wID[DEG_Limma_wID$logFC >0,]
Limma_neg <- DEG_Limma_wID[DEG_Limma_wID$logFC<0,]

saveRDS(Limma_pos , file = "Upreg_Limma.RDS")
saveRDS(Limma_neg, file = "Downreg_Limma.RDS" )

edgeR_pos <- DEG_edgeR_wID[DEG_edgeR_wID$logFC >0,]
edgeR_neg <- DEG_edgeR_wID[DEG_edgeR_wID$logFC<0,]

saveRDS(edgeR_pos , file = "Upreg_edgeR.RDS")
saveRDS(edgeR_neg, file = "Downreg_edgeR.RDS" )

#Creating venn diagrams between edgeR and limma

venndiagram <- venn.diagram(
  x = list(Limma_pos$ID, Limma_neg$ID,edgeR_pos$ID ,edgeR_neg$ID ),
  category.names = c("Upreg_Limma" , "Downreg_Limma","Upreg_edgeR" , "Downreg_edgeR"),
  filename = NULL
  
)
png('venndiagram.png', width = 20,height =14,res = 300, units = "cm")
grid.draw(venndiagram)
dev.off()

##########################################################################################################
#top 200 genes for up- and down regulation

Limma_pos200 <-Limma_pos[order(Limma_pos$logFC,decreasing = TRUE), ] [1:200,]
Limma_neg200 <- Limma_neg[order(Limma_neg$logFC), ] [1:200,]

edgeR_pos200 <-edgeR_pos[order(edgeR_pos$logFC,decreasing = TRUE), ] [1:200,]
edgeR_neg200 <- edgeR_neg[order(edgeR_neg$logFC), ] [1:200,]

#Create correlation plot 
png('correlationplot.png', width = 16,height =8,res = 300, units = "cm")
par(mfrow=c(1,2))
plot(Limma_pos200$logFC, edgeR_pos200$logFC ,xlab = "Limma logFC",ylab = "edgeR logFC",main = "Top 200 upregulated DEGs")
plot(Limma_neg200$logFC, edgeR_neg200$logFC ,xlab = "Limma logFC",ylab = "edgeR logFC",main = "Top 200 downregulated DEGs ")
dev.off()


##########################################################################################################
#Enrichment analysis

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021","KEGG_2021_Human",
         "Reactome_2016","MSigDB_Hallmark_2020")

#Enrichment of genes found by both methods
DEGs_pos <- intersect(Limma_pos$ID ,edgeR_pos$ID )
DEGs_neg <- intersect(Limma_neg$ID ,edgeR_neg$ID )

#Perform enrichment analysis
enriched_pos <- enrichr(DEGs_pos , dbs)
save(enriched_pos,file = "enriched_pos.rda")

enriched_neg <- enrichr(DEGs_neg , dbs)
save(enriched_neg,file = "enriched_neg.rda")

#Enrichment plot DEGs upregulated

png('Enrichment_MF_pos.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_pos[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_CC_pos.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_pos[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_BP_pos.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_pos[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_KEGG_pos.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_pos[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Re_pos.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_pos[[5]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Hallmark_pos.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_pos[[6]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#Enrichment plot DEGs downregulated
png('Enrichment_MF_neg.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_neg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_CC_neg.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_neg[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_BP_neg.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_neg[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_KEGG_neg.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_neg[[4]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Re_neg.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_neg[[5]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Hallmark_neg.png', width = 20,height =10,res = 350, units = "cm")
par(mfrow=c(1,2))
plotEnrich(enriched_neg[[6]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

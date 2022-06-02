library(TCGAbiolinks)
library(ggplot2)
library(ggridges)


#load files and create ID column for edgeR DEGs
DEG_Limma_wID <- get(load(file = "../DEG_Limma.rda" ))

DEG_edgeR_wID <- get(load(file = "../DEG_edgeR.rda" ))
DEG_edgeR_wID$ID <- rownames(DEG_edgeR_wID) 

#Find DEG IDs present in both methods
DEGs_Limmaboth <- subset(DEG_Limma_wID,ID %in% intersect(DEG_Limma_wID$ID ,DEG_edgeR_wID$ID))

DEGs_edgeRboth <- subset(DEG_edgeR_wID,ID %in% intersect(DEG_Limma_wID$ID ,DEG_edgeR_wID$ID))

##########################################################################################################################
#create volcano plots

#volcano for limma
label_limmaboth <- subset(DEGs_Limmaboth, -log10(adj.P.Val) > 50)
label_limmaboth <- label_limmaboth[order(label_limmaboth$logFC ),] 
maxDown_limmaboth <- label_limmaboth$ID [1:2] 
maxUp_limmaboth <- label_limmaboth$ID [(length(label_limmaboth$ID )-1):length(label_limmaboth$ID )] 

volcano_Limmaboth <- TCGAVisualize_volcano(x=DEGs_Limmaboth$logFC , y=DEGs_Limmaboth$adj.P.Val , names = DEGs_Limmaboth$ID,  
                                           show.names = "highlighted", highlight = c(maxDown_limmaboth, maxUp_limmaboth),
                                          highlight.color = "orange",  ylab = "-log10(P-value adjusted)", xlab = "logFC",
                                           title = "Volcano Limmaboth", filename = "volcano_Limmaboth.png" )

#volcano for edgeR
label_edgeRboth <- subset(DEGs_edgeRboth, -log10(FDR) > 50)
label_edgeRboth <- label_edgeRboth[order(label_edgeRboth$logFC ),] 
maxDown_edgeRboth <- label_edgeRboth$ID [1:2] 
maxUp_edgeRboth <- label_edgeRboth$ID [(length(label_edgeRboth$ID )-1):length(label_edgeRboth$ID )] 

volcano_edgeRboth <- TCGAVisualize_volcano(x=DEGs_edgeRboth$logFC, y=DEGs_edgeRboth$FDR, names = DEGs_edgeRboth$ID,  
                                           show.names = "highlighted", highlight = c(maxDown_edgeRboth, maxUp_edgeRboth),
                                           highlight.color = "orange", ylab = "-log10(FDR P-value adjusted)",xlab = "logFC",
                                           title = "Volcano edgeRboth", filename = "volcano_edgeRboth.png" )

######################################################################################
#Function for generation dataframe to be used for ridgeline plot
DEG_enrich_collection <- function(enrich_output, process_num, Limma_DEG ){
  gen_enrich <- strsplit(enrich_output$Genes[1:process_num], split = ";")
  total_genes <- unlist(gen_enrich)
  
  #Initate Terms as empty vector
  Terms <- c()
  for (Term in 1:process_num){
    Terms <- append(Terms,rep(substr(enrich_output$Term[Term],1,30),length(gen_enrich[[Term]] ) ))
    
  } 
  df_enrich <- data.frame(Term = Terms,
                          Genes = total_genes,
                          logFC = Limma_DEG$logFC[ match(total_genes,Limma_DEG$ID )],
                          Significance = -log10(Limma_DEG$adj.P.Val[ match(total_genes,Limma_DEG$ID )])
  )
  return(df_enrich)
  
} 
######################################################################################
#load information about enrichment analysis
enriched_pos <- get(load("../enriched_pos.rda"))
enriched_neg <- get(load("../enriched_neg.rda"))

enriched_posBP <- enriched_pos[["GO_Biological_Process_2021"]]
enriched_negBP <- enriched_neg[["GO_Biological_Process_2021"]]
enriched_posKEGG <- enriched_pos[["KEGG_2021_Human"]]
enriched_negKEGG <- enriched_neg[["KEGG_2021_Human"]]
enriched_posRE <- enriched_pos[["Reactome_2016"]]
enriched_negRE <- enriched_neg[["Reactome_2016"]]

BP_pos <- DEG_enrich_collection(enriched_posBP, 5,DEGs_Limmaboth )
BP_neg <- DEG_enrich_collection(enriched_negBP, 5,DEGs_Limmaboth )

KEGG_pos <- DEG_enrich_collection(enriched_posKEGG, 5,DEGs_Limmaboth )
KEGG_neg <- DEG_enrich_collection(enriched_negKEGG, 5,DEGs_Limmaboth )

RE_pos <- DEG_enrich_collection(enriched_posRE, 5,DEGs_Limmaboth )
RE_neg <- DEG_enrich_collection(enriched_negRE, 5,DEGs_Limmaboth )


df_enrich_merged <- do.call("rbind", list(RE_neg, KEGG_neg, BP_neg,RE_pos, KEGG_pos, BP_pos))

#Ridgeline plot with logFC
ggplot(df_enrich_merged, aes(x = logFC, y = Term, fill =stat(x))) +
  geom_density_ridges_gradient() + scale_fill_viridis_c(name = "logFC", option = "C")+
  theme(axis.text = element_text(size = 12),axis.title=element_text(size=14),legend.title=element_text(size=16),
        legend.text = element_text(size = 12))
ggsave("Ridgeline_logFC.png",width = 12, height = 10 )

#Ridgeline plot with -log10(adj.P.val)
ggplot(df_enrich_merged, aes(x = Significance, y = Term, fill =stat(x))) +
  geom_density_ridges_gradient() + scale_fill_viridis_c(name = "-log10(adj.P.Value)", option = "C") +
  scale_x_continuous(name = "-log10(adj.P.Value)")+
  theme(axis.text = element_text(size = 12),axis.title=element_text(size=14),legend.title=element_text(size=16),
        legend.text = element_text(size = 12))
ggsave("Ridgeline_adj.P.val.png",width = 12, height = 10 )


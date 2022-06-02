
#load libraries
library(grid)
library(gridBase)
library(enrichR)
library(ComplexHeatmap)

DEG_Limma_CN_LOW <- get(load(file = "DEG_CN_LOW_Limmaboth.rda" ))
DEG_Limma_CN_HIGH <- get(load(file = "DEG_CN_HIGH_Limmaboth.rda" ))
DEG_Limma_MSI <- get(load(file = "DEG_MSI_Limmaboth.rda" ))
DEG_Limma_POLE <- get(load(file = "DEG_POLE_Limmaboth.rda" ))

#Find up and down regulated
CN_LOW_pos <- DEG_Limma_CN_LOW[DEG_Limma_CN_LOW$logFC>0,]
CN_LOW_neg <- DEG_Limma_CN_LOW[DEG_Limma_CN_LOW$logFC<0,]

CN_HIGH_pos <- DEG_Limma_CN_HIGH[DEG_Limma_CN_HIGH$logFC>0,]
CN_HIGH_neg <- DEG_Limma_CN_HIGH[DEG_Limma_CN_HIGH$logFC<0,]

MSI_pos <- DEG_Limma_MSI[DEG_Limma_MSI$logFC>0,]
MSI_neg <- DEG_Limma_MSI[DEG_Limma_MSI$logFC<0,]

POLE_pos <- DEG_Limma_POLE[DEG_Limma_POLE$logFC>0,]
POLE_neg <- DEG_Limma_POLE[DEG_Limma_POLE$logFC<0,]


#Enrichment analysis

dbs <- c("GO_Biological_Process_2021")
enriched_Limma_CN_LOW_pos <- enrichr(CN_LOW_pos$ID , dbs)
enriched_Limma_CN_LOW_neg <- enrichr(CN_LOW_neg$ID , dbs)

enriched_Limma_CN_HIGH_pos <- enrichr(CN_HIGH_pos$ID , dbs)
enriched_Limma_CN_HIGH_neg <- enrichr(CN_HIGH_neg$ID , dbs)

enriched_Limma_MSI_pos <- enrichr(MSI_pos$ID , dbs)
enriched_Limma_MSI_neg <- enrichr(MSI_neg$ID , dbs)

enriched_Limma_POLE_pos <- enrichr(POLE_pos$ID , dbs)
enriched_Limma_POLE_neg <- enrichr(POLE_neg$ID , dbs)


#view results CN_LOW
png('Enrichment_Limma_CN_LOW_pos.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_CN_LOW_pos[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Limma_CN_LOW_neg.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_CN_LOW_neg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#view results CN_HIGH
png('Enrichment_Limma_CN_HIGH_pos.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_CN_HIGH_pos[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Limma_CN_HIGH_neg.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_CN_HIGH_neg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#view results MSI
png('Enrichment_Limma_MSI_pos.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_MSI_pos[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Limma_MSI_neg.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_MSI_neg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

#view results POLE
png('Enrichment_Limma_POLE_pos.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_POLE_pos[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()

png('Enrichment_Limma_POLE_neg.png', width = 20,height =10,res = 350, units = "cm")
plotEnrich(enriched_Limma_POLE_neg[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
dev.off()


###############################################################################################################
#Create Upset plots

#List with DEGs in different subtypes
listInput <- list("CN HIGH" =DEG_Limma_CN_HIGH$ID, "CN LOW" = DEG_Limma_CN_LOW$ID,
                  MSI = DEG_Limma_MSI$ID, POLE = DEG_Limma_POLE$ID  )

comb_mat <- make_comb_mat(listInput)

subtype <- c("CN HIGH","CN LOW","MSI","POLE")
png("Upset_plot.png", width = 16,height =21,res = 350, units = "cm")
UpSet(comb_mat, set_order = subtype, pt_size = unit(3, "mm"), lwd = 2, height = unit(4, "cm"),
      comb_col = c("orange","yellow","red","black")[comb_degree(comb_mat)],
      top_annotation = upset_top_annotation(comb_mat, 
                                            height = unit(15, "cm"),
                                            ylim = c(0,4000),
                                            bar_width = 0.7, 
                                            axis_param = list(side = "left", at = seq(0,4000,400)),
                                            annotation_name_side = "left", 
                                            annotation_name_gp = gpar(cex = 1), 
                                            annotation_name_offset = unit(1.5,"cm")),
      right_annotation = upset_right_annotation(comb_mat, 
                                                width = unit(4, "cm"), 
                                                gp = gpar(fill = "lightgreen"),
                                                axis_param = list(at = seq(0,6000,1000)), 
                                                annotation_name_offset = unit(1, "cm")),
      row_names_gp = gpar(fontsize = 12)
      )
dev.off()

#Create list with up and downregulated DEGs in subtypes
listInput_updown <- list("CN HIGH UP" =CN_HIGH_pos$ID,  "CN LOW UP" =CN_LOW_pos$ID,
                         "MSI UP" = MSI_pos$ID, "POLE UP" =POLE_pos$ID,
                         "CN HIGH DOWN" =CN_HIGH_neg$ID,  "CN LOW DOWN" =CN_LOW_neg$ID,
                         "MSI DOWN" = MSI_neg$ID, "POLE DOWN" =POLE_neg$ID)

comb_mat_updown <- make_comb_mat(listInput_updown)

subtype_updown <- c("CN HIGH UP","CN LOW UP","MSI UP","POLE UP","CN HIGH DOWN","CN LOW DOWN",
             "MSI DOWN","POLE DOWN")

png("Upset_UpDown_plot.png", width = 18,height =21,res = 350, units = "cm")
UpSet(comb_mat_updown,set_order = subtype_updown, pt_size = unit(3, "mm"), lwd = 3, height = unit(4, "cm"),
      comb_col = c( "coral","orange","yellow","red","darkred","blueviolet","darkolivegreen4","black")[comb_degree(comb_mat_updown)],
      top_annotation = upset_top_annotation(comb_mat_updown, 
                                            height = unit(15, "cm"),
                                            ylim = c(0,2000),
                                            bar_width = 0.7, 
                                            axis_param = list(side = "left", at = seq(0,2000,200)),
                                            annotation_name_side = "left", 
                                            annotation_name_gp = gpar(cex = 1), 
                                            annotation_name_offset = unit(1.5,"cm")),
      right_annotation = upset_right_annotation(comb_mat_updown, 
                                                width = unit(4, "cm"), 
                                                gp = gpar(fill = "lightgreen"),
                                                axis_param = list(at = seq(0,3500,1000)), 
                                                annotation_name_offset = unit(1, "cm")),
      row_names_gp = gpar(fontsize = 12)
)
dev.off()

#######################################################################################################
#Investigate gene duplicates
duplicate_CN_LOW <- intersect(CN_LOW_neg$ID, CN_LOW_pos$ID)
duplicate_CN_HIGH <- intersect(CN_HIGH_neg$ID, CN_HIGH_pos$ID)
duplicate_MSI <- intersect(MSI_neg$ID, MSI_pos$ID)
duplicate_POLE <- intersect(POLE_neg$ID, POLE_pos$ID)

#######################################################################################################
#Find genes that are the most Up- and downregulated

CN_HIGH_up <- subset(CN_HIGH_pos[which.max(CN_HIGH_pos$logFC),] , select = ID)
CN_LOW_up <- subset(CN_LOW_pos[which.max(CN_LOW_pos$logFC),] , select = ID)
MSI_up <- subset(MSI_pos[which.max(MSI_pos$logFC),] , select = ID)
POLE_up <- subset(POLE_pos[which.max(POLE_pos$logFC),] , select = ID)

CN_HIGH_down <- subset(CN_HIGH_neg[which.min(CN_HIGH_neg$logFC),] , select = ID)
CN_LOW_down <- subset(CN_LOW_neg[which.min(CN_LOW_neg$logFC),] , select = ID)
MSI_down <- subset(MSI_neg[which.min(MSI_neg$logFC),] , select = ID)
POLE_down <- subset(POLE_neg[which.min(POLE_neg$logFC),] , select = ID)

check_subtypes <- function(subtype_gene, subtype){
  check_CN_HIGH_pos<- ""
  check_CN_HIGH_neg<- ""
  check_CN_LOW_pos<- ""
  check_CN_LOW_neg<- ""
  check_MSI_pos<- ""
  check_MSI_neg<- ""
  check_POLE_pos<- ""
  check_POLE_neg <- ""
  
  if (subtype_gene %in% CN_HIGH_pos$ID ){
    check_CN_HIGH_pos <- "CN_HIGH_pos"
  }  
  if (subtype_gene %in% CN_HIGH_neg$ID ){
    check_CN_HIGH_neg <- "CN_HIGH_neg"
  }
  if (subtype_gene %in% CN_LOW_pos$ID ){
    check_CN_LOW_pos <- "CN_LOW_pos"
  }
  if (subtype_gene %in% CN_LOW_neg$ID ){
    check_CN_LOW_neg <- "CN_LOW_neg"
  }
  if (subtype_gene %in% MSI_pos$ID ){
    check_MSI_pos <- "MSI_pos"
  }
  if (subtype_gene %in% MSI_neg$ID ){
    check_MSI_neg <- "MSI_neg"
  }
  if (subtype_gene %in% POLE_pos$ID ){
    check_POLE_pos <- "POLE_pos"
  }
  if (subtype_gene %in% POLE_neg$ID ){
    check_POLE_neg <- "POLE_neg"
  }
  
  return(paste(subtype,"are found in",check_CN_HIGH_pos,check_CN_HIGH_neg,
               check_CN_LOW_pos,check_CN_LOW_neg,check_MSI_pos,
               check_MSI_neg,check_POLE_pos,check_POLE_neg))
  
} 

check_subtypes (CN_HIGH_up, "CN_HIGH_up")
check_subtypes (CN_HIGH_down, "CN_HIGH_down")
check_subtypes (CN_LOW_up, "CN_LOW_up")
check_subtypes (CN_LOW_down, "CN_LOW_down")
check_subtypes (MSI_up, "MSI_up")
check_subtypes (MSI_down, "MSI_down")
check_subtypes (POLE_up, "POLE_up")
check_subtypes (POLE_down, "POLE_down")

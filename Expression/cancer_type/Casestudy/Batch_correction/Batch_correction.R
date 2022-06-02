library(ggplot2)
library(TCGAbiolinks)
library(limma)

UCEC_dataFilt_HUGO <- get(load("../UCEC_dataFilt_HUGO_wID.rda"))
CN_LOW <- get(load("../Subtype/CN_LOW.rda"))
CN_HIGH <- get(load("../Subtype/CN_HIGH.rda"))
MSI <- get(load("../Subtype/MSI.rda"))
POLE <- get(load("../Subtype/POLE.rda"))

#First do log transformation with voom function
UCEC_dataFilt_HUGO_transform <- voom(UCEC_dataFilt_HUGO)
UCEC_dataFilt_HUGO_df <- as.data.frame(UCEC_dataFilt_HUGO_transform[["E"]])


#The MDSPlot function takes following arguments; Genecount matrix for all,
#genecount matrix for each subtype. It returns a PCA plot for all the samples
#with the subtypes and normal marked.

MDSPlot <- function(my.data, CN_LOW, CN_HIGH, MSI, POLE) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  
  #Do label with subtypes by TCGA barcode
  for (cell in 1:length(res$names)) {
    if (res$names[cell] %in% colnames(CN_LOW) ){
      res$Subtype[cell]  = "CN_LOW"
    }else if (res$names[cell] %in% colnames(CN_HIGH)){
      res$Subtype[cell]  = "CN_HIGH"
    }else if (res$names[cell] %in% colnames(MSI)){
      res$Subtype[cell]  = "MSI"
    }else if (res$names[cell] %in% colnames(POLE)){
      res$Subtype[cell]  = "POLE"
    }
    else res$Subtype[cell]  = "NORMAL"
  }     
  
  ggplot(res, aes(x=M1, y=M2, color = Subtype)) + geom_point(size = 3)+
    ggtitle("PCA UCEC Log transformed raw data") +
    
    theme_light()+
    theme(legend.text = element_text(size = 16), axis.title=element_text(size=16),legend.title=element_text(size=18),
          plot.title = element_text(size=20),axis.text=element_text(size=14))
  
}

#PCA plot not batch corrected
MDSPlot(UCEC_dataFilt_HUGO_df, CN_LOW,CN_HIGH, MSI, POLE)
ggsave(filename = "UCEC_non_batch_corrected.png")

UCEC_batch_corrected_Plate <- TCGAbatch_Correction(tabDF = UCEC_dataFilt_HUGO_transform$E , batch.factor="Plate", adjustment=NULL)

#PCA plot of Plate batch corrected data
MDSPlot(UCEC_batch_corrected_Plate, CN_LOW,CN_HIGH, MSI, POLE)
ggsave(filename = "UCEC_plate_batch_corrected.png")

#Error: need finite ylim values
UCEC_batch_corrected_TSS <- TCGAbatch_Correction(tabDF = UCEC_dataFilt_HUGO_transform$E , batch.factor="TSS", adjustment=NULL)

UCEC_batch_corrected_Portion <- TCGAbatch_Correction(tabDF = UCEC_dataFilt_HUGO_transform$E , batch.factor="Portion", adjustment=NULL)

# Error: contrasts can be applied only to factors with 2 or more levels
UCEC_batch_corrected_Center <- TCGAbatch_Correction(tabDF = UCEC_dataFilt_HUGO_transform$E , batch.factor="Sequencing Center", adjustment=NULL)

#Error: year data not provided
UCEC_batch_corrected_Year <- TCGAbatch_Correction(tabDF = UCEC_dataFilt_HUGO_transform$E , batch.factor="Year", adjustment=NULL)

#Explore the clinical data
my_IDs <- get_IDs(UCEC_dataFilt_HUGO_transform$E)
Summary_batch <- summary(my_IDs,  maxsum = 31)

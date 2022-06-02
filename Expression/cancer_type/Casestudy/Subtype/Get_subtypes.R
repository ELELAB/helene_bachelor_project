
library(TCGAbiolinks)
library(TCGAutils)

get_GDC_data <- function(cancer, subtype){
  
  # Find subtype specific barcodes
  pancancer_table <- PanCancerAtlas_subtypes()
  pancancer_cancer <- pancancer_table[which(pancancer_table$cancer.type == cancer),]
  barcodes <- pancancer_cancer$pan.samplesID[which(pancancer_cancer$Subtype_Selected == subtype)]
  
  
  # remove normal samples
  if(length(strsplit(barcodes,split="-")[[1]]) > 3){
    tp_barcodes <- TCGAquery_SampleTypes(barcodes, "TP")
    barcodes_three <- lapply(strsplit(tp_barcodes,split="-"), "[", c(1:3))
    subtype_barcodes <- sapply(barcodes_three, paste0, collapse="-")
  }else{
    
    subtype_barcodes <- barcodes} 
}

UCEC_dataFilt_HUGO <- get(load("../UCEC_dataFilt_HUGO_wID.rda"))

UCEC_submitterID <- TCGAbiospec(colnames(UCEC_dataFilt_HUGO))$submitter_id

#Find CN_LOW subtype
barcode_CN_LOW <- get_GDC_data ("UCEC","UCEC.CN_LOW")
match_CN_LOW <- na.omit(match(barcode_CN_LOW, UCEC_submitterID))

CN_LOW <- UCEC_dataFilt_HUGO[,match_CN_LOW] 
save(CN_LOW, file = "CN_LOW.rda")

#Find CN_HIGH subtype
barcode_CN_HIGH <- get_GDC_data ("UCEC","UCEC.CN_HIGH")
match_CN_HIGH <- na.omit(match(barcode_CN_HIGH, UCEC_submitterID))

CN_HIGH <- UCEC_dataFilt_HUGO[,match_CN_HIGH] 
save(CN_HIGH, file = "CN_HIGH.rda")

#Find MSI subtype
barcode_MSI <- get_GDC_data ("UCEC","UCEC.MSI")
match_MSI <- na.omit(match(barcode_MSI, UCEC_submitterID))

MSI <- UCEC_dataFilt_HUGO[,match_MSI] 
save(MSI, file = "MSI.rda")

#Find POLE subtype
barcode_POLE <- get_GDC_data ("UCEC","UCEC.POLE")
match_POLE <- na.omit(match(barcode_POLE, UCEC_submitterID))

POLE <- UCEC_dataFilt_HUGO[,match_POLE] 
save(POLE, file = "POLE.rda")

#dataframe with number of samples in each subtype
subtype_distribution <- data.frame(
                                  "CN_LOW"=length(match_CN_LOW),
                                  "CN_HIGH"=length(match_CN_HIGH),
                                  "MSI"=length(match_MSI),
                                  "POLE"=length(match_POLE),
                                  row.names = "Number of samples"
  
)



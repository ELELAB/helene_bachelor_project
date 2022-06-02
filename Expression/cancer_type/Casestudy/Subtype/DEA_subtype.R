
#load libaries
library(TCGAbiolinks)
library(limma)
library(edgeR)

#import revised DEA function
source("../../TCGAanalyze_DEA_revised.R")

#Load files for analysis
CN_LOW <- get(load(file = "CN_LOW.rda"))
CN_HIGH <- get(load(file = "CN_HIGH.rda"))
MSI <- get(load(file = "MSI.rda"))
POLE <- get(load(file = "POLE.rda"))

UCEC_dataFilt_HUGO_wID <- get(load("../UCEC_dataFilt_HUGO_wID.rda"))

########################################################################
#Analysis for HUGO
# selection of normal samples "NT" from HUGO
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(UCEC_dataFilt_HUGO_wID),
                                   typesample = c("NT"))


#DEA with Limma by revised DEA function of CN_LOW
dataDEGs_Limma_CN_LOW <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                          mat2 = CN_LOW,
                                          Cond1type = "Normal",
                                          Cond2type = "Tumor",
                                          fdr.cut = 0.1,
                                          logFC.cut = 1,
                                          pipeline = "limma",
                                          voom = TRUE,
                                          batch.factors ="TSS"
                                          
)
#DEA with edgeR GLM method by revised DEA function of CN_LOW
dataDEGs_edgeR_CN_LOW <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                                 mat2 = CN_LOW,
                                          Cond1type = "Normal",
                                          Cond2type = "Tumor",
                                          fdr.cut = 0.1,
                                          logFC.cut = 1,
                                          pipeline = "edgeR",
                                          method = "glmLRT",
                                          voom = FALSE,
                                          batch.factors ="TSS"
                                          
)
dataDEGs_edgeR_CN_LOW$ID <- rownames(dataDEGs_edgeR_CN_LOW)

CN_LOW_Limmaboth <- subset(dataDEGs_Limma_CN_LOW,ID %in% intersect(dataDEGs_Limma_CN_LOW$ID ,dataDEGs_edgeR_CN_LOW$ID ))

CN_LOW_edgeRboth <- subset(dataDEGs_edgeR_CN_LOW,ID %in% intersect(dataDEGs_Limma_CN_LOW$ID ,dataDEGs_edgeR_CN_LOW$ID ))

save(CN_LOW_Limmaboth, file = "DEG_CN_LOW_Limmaboth.rda" )
save(CN_LOW_edgeRboth, file = "DEG_CN_LOW_edgeRboth.rda" )

#DEA with Limma by revised DEA function of CN_HIGH
dataDEGs_Limma_CN_HIGH <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                                 mat2 = CN_HIGH,
                                                 Cond1type = "Normal",
                                                 Cond2type = "Tumor",
                                                 fdr.cut = 0.1,
                                                 logFC.cut = 1,
                                                 pipeline = "limma",
                                                 voom = TRUE,
                                                 batch.factors ="TSS"
                                                 
)
#DEA with edgeR GLM method by revised DEA function of CN_HIGH
dataDEGs_edgeR_CN_HIGH <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                                 mat2 = CN_HIGH,
                                                 Cond1type = "Normal",
                                                 Cond2type = "Tumor",
                                                 fdr.cut = 0.1,
                                                 logFC.cut = 1,
                                                 pipeline = "edgeR",
                                                 method = "glmLRT",
                                                 voom = FALSE,
                                                 batch.factors ="TSS"
                                                 
)
dataDEGs_edgeR_CN_HIGH$ID <- rownames(dataDEGs_edgeR_CN_HIGH)

CN_HIGH_Limmaboth <- subset(dataDEGs_Limma_CN_HIGH,ID %in% intersect(dataDEGs_Limma_CN_HIGH$ID ,dataDEGs_edgeR_CN_HIGH$ID ))

CN_HIGH_edgeRboth <- subset(dataDEGs_edgeR_CN_HIGH,ID %in% intersect(dataDEGs_Limma_CN_HIGH$ID ,dataDEGs_edgeR_CN_HIGH$ID ))

save(CN_HIGH_Limmaboth, file = "DEG_CN_HIGH_Limmaboth.rda" )
save(CN_HIGH_edgeRboth, file = "DEG_CN_HIGH_edgeRboth.rda" )


#DEA with Limma by revised DEA function of MSI
dataDEGs_Limma_MSI <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                                 mat2 = MSI,
                                                 Cond1type = "Normal",
                                                 Cond2type = "Tumor",
                                                 fdr.cut = 0.1,
                                                 logFC.cut = 1,
                                                 pipeline = "limma",
                                                 voom = TRUE,
                                                 batch.factors ="TSS"
                                                 
)
#DEA with edgeR GLM method by revised DEA function of MSI
dataDEGs_edgeR_MSI <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                                  mat2 = MSI,
                                                  Cond1type = "Normal",
                                                  Cond2type = "Tumor",
                                                  fdr.cut = 0.1,
                                                  logFC.cut = 1,
                                                  pipeline = "edgeR",
                                                  method = "glmLRT",
                                                  voom = FALSE,
                                                  batch.factors ="TSS"
                                                  
)
dataDEGs_edgeR_MSI$ID <- rownames(dataDEGs_edgeR_MSI)

MSI_Limmaboth <- subset(dataDEGs_Limma_MSI,ID %in% intersect(dataDEGs_Limma_MSI$ID ,dataDEGs_edgeR_MSI$ID ))

MSI_edgeRboth <- subset(dataDEGs_edgeR_MSI,ID %in% intersect(dataDEGs_Limma_MSI$ID ,dataDEGs_edgeR_MSI$ID ))

save(MSI_Limmaboth, file = "DEG_MSI_Limmaboth.rda" )
save(MSI_edgeRboth, file = "DEG_MSI_edgeRboth.rda" )

#DEA with Limma by revised DEA function of POLE
dataDEGs_Limma_POLE <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                                 mat2 = POLE,
                                                 Cond1type = "Normal",
                                                 Cond2type = "Tumor",
                                                 fdr.cut = 0.1,
                                                 logFC.cut = 1,
                                                 pipeline = "limma",
                                                 voom = TRUE,
                                                 batch.factors ="TSS"
                                                 
)
#DEA with edgeR GLM method by revised DEA function of MSI
dataDEGs_edgeR_POLE <- TCGAanalyze_DEA_revised(mat1 = UCEC_dataFilt_HUGO_wID[,samplesNT] ,
                                              mat2 = POLE,
                                              Cond1type = "Normal",
                                              Cond2type = "Tumor",
                                              fdr.cut = 0.1,
                                              logFC.cut = 1,
                                              pipeline = "edgeR",
                                              method = "glmLRT",
                                              voom = FALSE,
                                              batch.factors ="TSS"
                                              
)
dataDEGs_edgeR_POLE$ID <- rownames(dataDEGs_edgeR_POLE)

POLE_Limmaboth <- subset(dataDEGs_Limma_POLE,ID %in% intersect(dataDEGs_Limma_POLE$ID ,dataDEGs_edgeR_POLE$ID ))

POLE_edgeRboth <- subset(dataDEGs_edgeR_POLE,ID %in% intersect(dataDEGs_Limma_POLE$ID ,dataDEGs_edgeR_POLE$ID ))

save(POLE_Limmaboth, file = "DEG_POLE_Limmaboth.rda" )
save(POLE_edgeRboth, file = "DEG_POLE_edgeRboth.rda" )



#This script performs a differential expression analysis and compare the
#differential expressed genes (DEGs) originating from two different normalization
#methods Limma-voom and edgeR GLM.

#load libaries
library(TCGAbiolinks)
library(limma)
library(edgeR)


# Get cancer as user input from command line
cancer <- as.character(commandArgs(trailingOnly = TRUE))

#import revised DEA function
source("../../TCGAanalyze_DEA_revised.R")

#Load files for analysis
dataFilt_HUGO <- get(load(file = paste0("../",cancer,"_dataFilt_HUGO.rda")))

#Filter for genes without ID
dataFilt_HUGO_wID <- dataFilt_HUGO[-which(rownames(dataFilt_HUGO) == ""),] 
save(dataFilt_HUGO_wID, file = paste0("../",cancer,"_dataFilt_HUGO_wID.rda"))

########################################################################
#Analysis for HUGO
# selection of normal samples "NT" from HUGO
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_HUGO_wID),
                                   typesample = c("NT"))

# selection of tumor samples "TP" from HUGO
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_HUGO_wID), 
                                   typesample = c("TP"))


#DEA with Limma by revised DEA function
dataDEGs_Limma <- TCGAanalyze_DEA_revised(mat1 = dataFilt_HUGO_wID[,samplesNT] ,
                                        mat2 = dataFilt_HUGO_wID[,samplesTP],
                                        Cond1type = "Normal",
                                        Cond2type = "Tumor",
                                        fdr.cut = 0.1,
                                        logFC.cut = 1,
                                        pipeline = "limma",
                                        voom = TRUE,
                                        batch.factors ="TSS"
                                        
)
save(dataDEGs_Limma, file = "DEG_Limma.rda" )


#DEA with edgeR GLM method by revised DEA function
dataDEGs_edgeR <- TCGAanalyze_DEA_revised(mat1 = dataFilt_HUGO_wID[,samplesNT],
                                          mat2 = dataFilt_HUGO_wID[,samplesTP],
                                          Cond1type = "Normal",
                                          Cond2type = "Tumor",
                                          fdr.cut = 0.1,
                                          logFC.cut = 1,
                                          pipeline = "edgeR",
                                          method = "glmLRT",
                                          voom = FALSE,
                                          batch.factors ="TSS"
                                          
)
save(dataDEGs_edgeR, file = "DEG_edgeR.rda" )


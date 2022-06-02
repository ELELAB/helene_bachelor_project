#This script performs a differential expression analysis and compare the
#differential expressed genes (DEGs) originating from two different normalization
#methods

#load libaries
library(TCGAbiolinks)
library(limma)

#import revised DEA function
source("../TCGAanalyze_DEA_revised.R")

#Load files for analysis
load(file = "../BRCA_getGLGC/BRCA_dataFilt.rda")
dataFilt_HUGO <- get(load(file = "../BRCA_getGLGC/BRCA_dataFilt_HUGO.rda"))
dataFilt_oldGCn <- get(load(file = "../BRCA_getGLGC/BRCA_dataFilt_oldGCn.rda"))
dataFilt_oldGCn_HUGO <- get(load(file = "../BRCA_getGLGC/BRCA_dataFilt_oldGCn_HUGO.rda"))

########################################################################
#Analysis for ENSEMBL
# selection of normal samples "NT" from ENSEMBL
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP" from ENSEMBL
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))


#DEA by revised DEA function
dataDEGs_ensembl <- TCGAanalyze_DEA_revised(mat1 = dataFilt[,samplesNT],
                                        mat2 = dataFilt[,samplesTP],
                                        Cond1type = "Normal",
                                        Cond2type = "Tumor",
                                        fdr.cut = 0.05,
                                        logFC.cut = 0.5,
                                        pipeline = "limma",
                                        voom = TRUE,
                                        batch.factors ="TSS"
                                        
)
save(dataDEGs_ensembl, file = "DEG_ensembl.rda" )

# selection of normal samples "NT" from oldGCn ENSEMBL 
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_oldGCn),
                                   typesample = c("NT"))

# selection of tumor samples "TP" from oldGCn ENSEMBL
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_oldGCn), 
                                   typesample = c("TP"))


#DEA by revised DEA function
dataDEGs_ensembl_oldGCn <- TCGAanalyze_DEA_revised(mat1 = dataFilt_oldGCn[,samplesNT],
                                        mat2 = dataFilt_oldGCn[,samplesTP],
                                        Cond1type = "Normal",
                                        Cond2type = "Tumor",
                                        fdr.cut = 0.05,
                                        logFC.cut = 0.5,
                                        pipeline = "limma",
                                        voom = TRUE,
                                        batch.factors ="TSS"
                                        
)
save(dataDEGs_ensembl_oldGCn, file = "DEG_ensembl_oldGCn.rda" )


########################################################################
#Analysis for HUGO
# selection of normal samples "NT" from HUGO
samplesNT_hugo <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_HUGO),
                                   typesample = c("NT"))

# selection of tumor samples "TP" from HUGO
samplesTP_hugo <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_HUGO), 
                                   typesample = c("TP"))


#DEA by revised DEA function
dataDEGs_hugo <- TCGAanalyze_DEA_revised(mat1 = dataFilt_HUGO[,samplesNT_hugo],
                                        mat2 = dataFilt_HUGO[,samplesTP_hugo],
                                        Cond1type = "Normal",
                                        Cond2type = "Tumor",
                                        fdr.cut = 0.05,
                                        logFC.cut = 0.5,
                                        pipeline = "limma",
                                        voom = TRUE,
                                        batch.factors ="TSS"
                                        
)
save(dataDEGs_hugo, file = "DEG_HUGO.rda" )

# selection of normal samples "NT" from oldGCn HUGO
samplesNT_oldGCn_hugo <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_oldGCn_HUGO),
                                        typesample = c("NT"))

# selection of tumor samples "TP" from oldGCN HUGO
samplesTP_oldGCn_hugo <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_oldGCn_HUGO), 
                                        typesample = c("TP"))


#DEA by revised DEA function
dataDEGs_oldGCn_hugo <- TCGAanalyze_DEA_revised(mat1 = dataFilt_oldGCn_HUGO[,samplesNT_oldGCn_hugo],
                                         mat2 = dataFilt_oldGCn_HUGO[,samplesTP_oldGCn_hugo],
                                         Cond1type = "Normal",
                                         Cond2type = "Tumor",
                                         fdr.cut = 0.05,
                                         logFC.cut = 0.5,
                                         pipeline = "limma",
                                         voom = TRUE,
                                         batch.factors ="TSS"
                                         
)
save(dataDEGs_oldGCn_hugo, file = "DEG_HUGO_oldGCn.rda" )

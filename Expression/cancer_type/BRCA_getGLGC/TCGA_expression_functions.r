#----------------- Functions TCGA gene expression data -----------------------------------------------

# Get TCGA data from GDC
#get_GDC_data <- function(cancer){
  
  #dir.create("data")
  #project <- paste0("TCGA-", cancer)

  # For BRCA, only include female samples
  #if(cancer == "BRCA"){
    #dataClin <- GDCquery_clinic(project = project, type = "clinical")
    #female_barcodes <- dataClin$submitter_id[which(dataClin$gender == "female")]

    #query.exp <- GDCquery(project = project,
                        #data.category = "Transcriptome Profiling",
                        #data.type = "Gene Expression Quantification",
                        #workflow.type = "HTSeq - Counts",
                        #sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                        #legacy = FALSE,
                        #barcode = female_barcodes)
  #}else{

    #query.exp <- GDCquery(project = project,
                        #data.category = "Transcriptome Profiling",
                        #data.type = "Gene Expression Quantification",
                        #workflow.type = "HTSeq - Counts",
                        #sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                        #legacy = FALSE)
  #}
  
  #GDCdownload(query.exp, method = "api")
  
  #tumor_exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = paste0("data/",cancer, ".exp.rda"), directory = "../../GDCdata")
  #return(tumor_exp)

#}

# Prepare samples
prepare_data <- function(gdc){
  dataPrep <- TCGAanalyze_Preprocessing(object = gdc,
                                        cor.cut = 0.6,
                                        filename = paste0("dataPrep_",cancer, "_array_array.png"))
  
  return(dataPrep)
} 

# Normalize samples
normalize_data_oldGCn <- function(dataPrep){
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneInfoHT,
                                        method = "gcContent")
 
  return(dataNorm)
}

# Normalize samples by new GC txt table and not the new rda table geneInfoHT, because 
#it will replace the auto-generated table used for normalize_data_oldGCn
geneinfo_GLGC <- read.table(file = "../../get_GLGC_content/GLGC_table.txt")
normalize_data <- function(dataPrep){
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneinfo_GLGC,
                                        method = "gcContent")
  
  return(dataNorm)
}



# Filter samples
filter_data <- function(dataNorm){
  dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                    method = "quantile",
                                    qnt.cut =  0.25)
  return(dataFilt)
}

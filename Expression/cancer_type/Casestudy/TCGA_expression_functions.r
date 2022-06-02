#----------------- Functions TCGA gene expression data -----------------------------------------------

# Prepare samples
prepare_data <- function(gdc){
  dataPrep <- TCGAanalyze_Preprocessing(object = gdc,
                                        cor.cut = 0.6,
                                        filename = paste0("dataPrep_",cancer, "_array_array.png"))
  
  return(dataPrep)
} 


# Normalize samples by new GC table
geneInfoHT <- get(load(file = "../../get_GLGC_content/geneInfoHT.rda"))
normalize_data <- function(dataPrep){
  dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                        geneInfo = geneInfoHT,
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

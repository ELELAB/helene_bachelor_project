
#------------------ TCGA data for mRNA expression ----------------------------


library(SummarizedExperiment)
library(TCGAbiolinks)
library(TCGAutils)
library(ggplot2)
library(ggrepel)
library(limma)

# Get cancer as user input from command line
cancer <- as.character(commandArgs(trailingOnly = TRUE))

# Source functions
source("TCGA_expression_functions.r")
source("TCGA_replicate_functions.r")

#read in data prepared by GDCprepare
gdc <- get(load(file = paste0("../../../data/",cancer,".exp.rda")))

#remove samples with read count sum less than 20 million reads
gdc_countmatrix <- assay(gdc)
samplesum <- colSums(gdc_countmatrix)

keep_samples_non_lowcount <- samplesum > 20*10**6
gdc_countmatrix_non_lowcount <- gdc_countmatrix[,keep_samples_non_lowcount] 

#################################################################################
#Explore technical replicates
#Running function for tumor samples TP
info_replicate_TP <- check_for_replicates(gdc_countmatrix_non_lowcount, cancer,"TP")

#Running function for normal samples NT
info_replicate_NT <- check_for_replicates(gdc_countmatrix_non_lowcount, cancer,"NT")

writeLines(c(paste("Logfile of replicates for",cancer,"\nTP = tumor, NT = normal\n "),info_replicate_TP,info_replicate_NT), paste0(cancer,"_logfile.txt"))

# Save as Ranged Summarized Experiment object for input to replicate analysis function
gdc_non_lowcount <- gdc[,keep_samples_non_lowcount] 

#Return samples to discard from data
discard_samples_TP <- replicate_analysis(gdc_non_lowcount, cancer, "TP")
discard_samples_NT <- replicate_analysis(gdc_non_lowcount, cancer, "NT")
discard_samples <- c(discard_samples_TP,discard_samples_NT)

keep_samples_non_replicate <-  !(colnames( gdc_non_lowcount) %in% discard_samples)

# Save new gdc without replicate as Ranged Summarized Experiment object for input to dataPrep function
gdc_non_replicate <- gdc_non_lowcount[,keep_samples_non_replicate] 

#################################################################################
# Pre process data
dataPrep <- prepare_data(gdc_non_replicate)
save(dataPrep, file = paste0(cancer, "_dataPrep.rda"))

# Save one additional copy, with HUGO names as row names instead of ensembl IDs
dataPrep_hugo <- dataPrep  
rownames(dataPrep_hugo) <- rowData(gdc)[match(rownames(dataPrep_hugo), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
save(dataPrep_hugo, file = paste0(cancer, "_dataPrep_HUGO.rda"))

# Normalize data (GC-count normalization)
dataNorm <- normalize_data(dataPrep)
save(dataNorm,file=paste0(cancer, "_dataNorm.rda"))  

# Save one additional copy, with HUGO names as row names instead of ensembl IDs
dataNorm_hugo <- dataNorm
rownames(dataNorm_hugo) <- rowData(gdc)[match(rownames(dataNorm_hugo), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
save(dataNorm_hugo, file = paste0(cancer, "_dataNorm_HUGO.rda"))

# Filter data
dataFilt <- filter_data(dataNorm)
save(dataFilt, file=paste0(cancer,"_dataFilt.rda"))

# Save one additional copy, with HUGO names as row names instead of ensembl IDs
dataFilt_hugo <- dataFilt
rownames(dataFilt_hugo) <- rowData(gdc)[match(rownames(dataFilt_hugo), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
save(dataFilt_hugo, file = paste0(cancer, "_dataFilt_HUGO.rda"))


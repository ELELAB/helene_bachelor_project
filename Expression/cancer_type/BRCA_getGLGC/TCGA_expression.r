
#------------------ TCGA data for mRNA expression ----------------------------

library(SummarizedExperiment)
library(TCGAbiolinks)

# Get cancer as user input from command line
cancer <- as.character(commandArgs(trailingOnly = TRUE))

# Source functions
source("TCGA_expression_functions.r")

# Original TCGA mRNA expression data was retrieved by this function but GDC update
#changed the parameters, so it was not comparable with the work already done
#gdc <- get_GDC_data(cancer)

#read in data prepared by GDCprepare
gdc <- get(load(file = paste0("../../../data/",cancer,".exp.rda")))

# Pre process data
dataPrep <- prepare_data(gdc)
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

# Normalize data (GC-count normalization with old table)
dataNorm_old <- normalize_data_oldGCn(dataPrep)
save(dataNorm_old,file=paste0(cancer, "_dataNorm_oldGCn.rda"))  

# Save one additional copy, with HUGO names as row names instead of ensembl IDs
dataNorm_old_hugo <- dataNorm_old
rownames(dataNorm_old_hugo) <- rowData(gdc)[match(rownames(dataNorm_old_hugo), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
save(dataNorm_old_hugo, file = paste0(cancer, "_dataNorm_oldGCn_HUGO.rda"))


# Filter data
dataFilt <- filter_data(dataNorm)
save(dataFilt, file=paste0(cancer,"_dataFilt.rda"))

# Save one additional copy, with HUGO names as row names instead of ensembl IDs
dataFilt_hugo <- dataFilt
rownames(dataFilt_hugo) <- rowData(gdc)[match(rownames(dataFilt_hugo), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
save(dataFilt_hugo, file = paste0(cancer, "_dataFilt_HUGO.rda"))

#Filter oldGCn data
dataFilt_old <- filter_data(dataNorm_old)
save(dataFilt_old, file=paste0(cancer,"_dataFilt_oldGCn.rda"))

# Save one additional copy, with HUGO names as row names instead of ensembl IDs
dataFilt_old_hugo <- dataFilt_old
rownames(dataFilt_old_hugo) <- rowData(gdc)[match(rownames(dataFilt_old_hugo), rowData(gdc)[,"ensembl_gene_id"]),"external_gene_name"]
save(dataFilt_old_hugo, file = paste0(cancer, "_dataFilt_oldGCn_HUGO.rda"))


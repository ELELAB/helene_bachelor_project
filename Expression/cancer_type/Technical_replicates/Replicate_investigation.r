
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
source("Replicate_functions.r")

#read in data prepared by GDCprepare
gdc <- get(load(file = paste0("../../../data/",cancer,".exp.rda")))
gdc_countmatrix <- assay(gdc)


#################################################################################
#Explore technical replicates
#Running function for tumor samples TP
info_replicate_TP <- check_for_replicates(gdc_countmatrix, cancer,"TP")

#Running function for normal samples NT
info_replicate_NT <- check_for_replicates(gdc_countmatrix, cancer,"NT")

writeLines(c(paste("Logfile of replicates for",cancer,"\nTP = tumor, NT = normal\n "),info_replicate_TP,info_replicate_NT), paste0(cancer,"_logfile.txt"))

#Return samples to discard from data
analysis_samples_TP <- replicate_analysis(gdc, cancer, "TP")
analysis_samples_NT <- replicate_analysis(gdc, cancer, "NT")


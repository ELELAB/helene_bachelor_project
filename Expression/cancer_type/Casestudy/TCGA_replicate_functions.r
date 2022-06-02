
#This script includes functions to extract replicates and analyze replicate data by generating
#barplots, coverage and PCA plots 

##########################################################################################################

#The check_for_replicates function takes a dataframe with gene count (samples as columns, gene as rows),
#cancer type and sample type and save a csv file containing TCGA_barcode, TCGA_submitterID and Match_id as column names.
#The function returns a character string explaining the number of unique and total replicates or if no replicates are found.

check_for_replicates <- function(tumor_exp, cancer, sample_type){

    #Get samples for sample_type
    samples <- TCGAquery_SampleTypes(barcode = colnames(tumor_exp), 
                                     typesample = sample_type)
  
    #Save submitter ID as patient variable
    patient <- TCGAbiospec(samples)$submitter_id 
    
    #save submitter ID of replicated 
    rep_patient <- patient[duplicated(patient) | duplicated(patient, fromLast=TRUE)]
    
    #Find unique replicates
    rep_patient_unique <-unique(rep_patient)
    
    #Create tables with duplicate information
    duplicate_table <- data.frame("TCGA_barcode" = samples,
                                  "TCGA_submitterID" = patient,
                                  "Match_id" = match(patient,rep_patient))
    write.csv(duplicate_table[order(duplicate_table$Match_id),], file = paste0(cancer,"_duplicate_info_",sample_type,".csv"))
    
    #Finding number of total and unique replicates
    Num_total_rep <- length(rep_patient)
    Num_unique_rep <- length(rep_patient_unique)
    
    if (Num_unique_rep == 0){
      replicate <- paste("No replicates are found in the",sample_type, "samples in this cancer type")
    } else{
      replicate <- paste("The number of unique replicates for the", sample_type, "samples are:",Num_unique_rep,
                         "\nThe Number of total replicates for the", sample_type, "samples are:", Num_total_rep)
    }
    return(replicate)
  
}


##########################################################################################################

#The genecount_subID functions saves the input gene count matrices with submitter ID (patient ID)
#as column names and return this matrix

get_genecount_subID <- function(genecount){
  data_subID <- genecount
  colnames(data_subID)<-TCGAbiospec(colnames(data_subID))$submitter_id
  
  return(data_subID)
} 

##########################################################################################################
#The get_match_id function takes following arguments which are for all the samples (unless stated);
#Gene count matrix, complete TCGA barcode, sample-plate identifier of barcode as character, all TCGA submitter IDs,
#and unique submitter ID for replicates.
#Returns a dataframe with Genecount, TCGA_barcode, Identifier, TCGA_submitter_ID, Samples, match_ID and Replicate as columns.

get_match_id <- function(genecount.all, barcode, barcode_label, all_submitterID, rep_patientID_unique){
  
  #calculate gene count sum for each sample and generate data frame 
  samplesum <- colSums(genecount.all)
  
  datas <- data.frame("Genecount" = samplesum,
                      "TCGA_barcode" = barcode,
                      "Identifier" = barcode_label,
                      "TCGA_submitter_ID" = all_submitterID)
  
  #Sort the sample sums numerically and match the submitter IDs to the replicate unique submitter_ids
  datas <- datas[order(datas$Genecount ),] 
  datas$Samples <-  seq(1:length(samplesum))
  datas$Match_id <- match(datas$TCGA_submitter_ID ,rep_patientID_unique)
  
  #Replace Non-replicates Match_id = NA with Match_id = 0
  datas[is.na(datas)] = 0
  for (cell in datas$Samples) {
    if (datas$Match_id[cell]  == 0 ){
      datas$Replicate[cell]  = "nonrep"
    }else{ 
      datas$Replicate[cell]  = "rep"
    }
  }
  
  return (datas)}


##########################################################################################################
#The get_match_id_rep function takes following arguments which are only for the replicates;
#Gene count matrix, TCGA barcode as character,submitter ID, unique submitter ID.
#Returns a dataframe with Genecount, TCGA_barcode and Match_id as columns.

get_match_id_rep <- function(rep_genecount, rep_barcode, rep_patientID, rep_patientID_unique){
  
  #calculate genecount sum for each sample and match the calculation to the unique submitter_ids
  samplesum <- colSums(rep_genecount)
  matches <-  match(names(samplesum),rep_patientID_unique)
  
  #Collect information in dataframe
  datas <- data.frame("Genecount" = samplesum,
                      "TCGA_barcode" = rep_barcode,
                      "Match_id" = matches)
  
  return (datas) 
}  

##########################################################################################################

#The barplotsum function takes following arguments which are only for the replicates;
#Dataframe consisting of Genecount, TCGA_barcode, Identifier, TCGA_submitter_ID,
#Samples, match_ID and Replicate as columns. And unique submitter_id, a color vector, cancer type and datatype.
#Returns a barplot of the replicates gene count sums, color coded by replicate.

barplotsum <- function(datas,rep_patientID_unique, my.cols, cancertype, datatype){
  
  #Generate barplot
  plots<- ggplot(data = datas ,aes(x= TCGA_barcode,y=Genecount,fill = as.factor(Match_id))) + 
    
    geom_bar(stat="identity") +
    ggtitle(paste("Coverage analysis for", cancertype, datatype))+
    scale_fill_manual(values=my.cols[2:(length(as.character(rep_patientID_unique))+1)], name = "Submitter ID",
                      labels = as.character(rep_patientID_unique)) +
    
    theme_light()+
    theme(axis.text.x = element_text(angle = 90,size=14),legend.text = element_text(size = 14),
          legend.title=element_text(size=16), plot.title = element_text(size=20),axis.title = element_text(size=16),
          axis.text.y = element_text(size=14),axis.text = element_text(color="black", face = "bold"))
  
  return(plots)
}  

##########################################################################################################

#The distribution_plot function takes following arguments which are for all the samples;
#Dataframe consisting of Genecount, TCGA_barcode, Identifier, TCGA_submitter_ID, Samples, match_ID
#and Replicate as columns. And unique submitter_id, a color vector, cancer type and datatype.
#It returns a scatterplot for all the samples gene count sums, with the replicates marked and color coded.

distribution_plot <- function(datas,rep_patientID_unique, my.cols, cancertype, datatype){
  
  #Generate distribution plot
  plots<- ggplot(data = datas ,aes(Samples,Genecount, color = as.factor(Match_id))) + 
    
    geom_point(aes(size = Replicate)) +
    ggtitle(paste("Coverage analysis for", cancertype, datatype )) +
    scale_colour_manual(values=my.cols[1:(length(as.character(rep_patientID_unique))+1)] , name = "Submitter ID",
                        labels = c("Non replicates",as.character(rep_patientID_unique))) + 
    
    theme_light()+
    theme(legend.text = element_text(size = 14),axis.text.x = element_blank(), axis.text.y = element_text(size = 14, color = "black"),
            axis.title=element_text(size=16), legend.title=element_text(size=16), plot.title = element_text(size=20)) + 
    
    geom_text_repel(data = subset(datas, Match_id>0), aes(Samples, Genecount, label = Identifier),vjust =-4,hjust=0.7, 
                    size = 5,color="black") +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  return(plots)
} 

###################################################################################################

#The MDSPlot function takes following arguments which are for all the samples;
#Genecount matrix, dataframe consisting of Genecount, TCGA_barcode, Identifier, TCGA_submitter_ID, Samples, match_ID
#and Replicate as columns. And Unique submitter_id, a color vector, cancer type and datatype
#It returns a PCA plot for all the samples with the replicates marked

MDSPlot <- function(my.data, match_table, my.labels, my.cols, cancertype, datatype) {
  my.data_transform <- voom(my.data)
  my.data_df <- as.data.frame(my.data_transform[["E"]])
  d<-dist(t(my.data_df))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  
  res$Label <- match_table$TCGA_submitter_ID[match(res$names,match_table$TCGA_barcode  )]  
  res$Identifier <-  match_table$Identifier[match(res$names,match_table$TCGA_barcode  )]  
  res$Match_id <- match_table$Match_id[match(res$names,match_table$TCGA_barcode  )]
  res$Replicate <- match_table$Replicate[match(res$names,match_table$TCGA_barcode  )]
  
  ggplot(res, aes(x=M1, y=M2, color = as.factor(Match_id))) + 
    
    geom_point(aes(size = Replicate))+ 
    ggtitle(paste("PCA", cancertype,"Log transformed", datatype))+
    scale_colour_manual(values=my.cols[1:(length(my.labels)+1)] , name = "Submitter ID",
                        labels = c("Non replicates",my.labels)) +
    theme_light()+
    
    theme(legend.text = element_text(size = 16), axis.title=element_text(size=16),legend.title=element_text(size=18),
          plot.title = element_text(size=20),axis.text=element_text(size=14)) +
    
    geom_text_repel(data = subset(res, Match_id>0),aes(x=M1,y=M2, label= Identifier),vjust =-5,hjust=1, 
                    size = 4.5,color="black") +
    guides(color = guide_legend(override.aes = list(size = 5))) 
  
}

#############################################################################################
#The replicate_analysis function takes a Ranged Summarized Experiment object, cancer type and sample type as input.
#The function finds for each unique replicate the sample with the highest gene count which we wish to keep and
#returns a list of samples to discard from the data.
#If no replicates are found an empty character vector is returned so no samples are removed and statement describing
#no replicates are found for this sample type are printed to the terminal.
#The function generate barplot and distribution plot with the functions barplotsum and distribution_plot

replicate_analysis <- function (gdc, cancer, sample_type) {
  #Extract genecount matrix
  genecount <- assay(gdc)
  
  #Extract samples for sample_type
  samples <- TCGAquery_SampleTypes(barcode = colnames(gdc), 
                                   typesample = sample_type)
  genecount <- genecount[,samples]
  save(genecount, file = paste0(cancer,"_genecountRaw_",sample_type,".rda"))
  
  #Extract information about barcodes
  barcode_info <- TCGAbiospec(samples)
  barcode_label <- substring(samples,14,25)
  all_submitterID <- barcode_info$submitter_id 
  
  #read in table with replicate information
  rep_info <- read.csv(file = paste0(cancer,"_duplicate_info_",sample_type,".csv"))
  if (all(is.na(rep_info$Match_id )) == FALSE){
    
    #find the replicates only by excluding NA data from rep_info_TP
    rep_patient_data <- na.omit(rep_info)
    rep_barcode <- as.character(rep_patient_data$TCGA_barcode)
    
    #find unique patient names
    rep_patient_uniqueID <- as.character(unique(rep_patient_data$TCGA_submitterID))
    
    #save gene count matrices with submitter ID as column names
    genecount_subID <- get_genecount_subID (genecount)
    
    #save only replicate gene count_matrix with submitter ID as column names
    rep_genecount_subID <- get_genecount_subID(genecount[,as.character(rep_patient_data$TCGA_barcode)] )
    
    #Generate match_id_table for all data and only for replicates
    match_id_table <- get_match_id(genecount_subID, samples, barcode_label, all_submitterID, rep_patient_uniqueID)
    
    match_id_table_rep <- get_match_id_rep(rep_genecount_subID,rep_barcode, rep_patient_data$TCGA_submitterID, rep_patient_uniqueID)
    
    #Find the highest gene count sum among each patient
    max_samples <- c()
    for (number in 1:max(match_id_table_rep$Match_id)){
      
      #Look through gene count sums for each patient
      table <- subset(match_id_table_rep, Match_id == number)
      
      #Save barcode for the sample with highest gene count sum
      max_sample <- table[which.max(table$Genecount),]
      max_samples[number]  <- as.character( max_sample$Sample_ID )
      
    }  
    
    #Save barcodes that are not max for their patient, to discard
    discard_samples <- setdiff(as.character(match_id_table_rep$Sample_ID), max_samples )
    
    #Define color vector for plots
    mycols = c("grey90", "coral2", "cadetblue2", "deeppink2","darkorchid1","khaki3", "yellow1",
               "salmon2", "orange1","olivedrab4", "plum", "grey22", "tomato3", "green2", "blue5")
    
    #Generate barplot with genecount sums for replicate samples
    barplotsum(match_id_table_rep, rep_patient_uniqueID, mycols, cancer, "raw data")
    ggsave(paste0(cancer,"_raw_barplot_",sample_type,".png"),width = 12, height = 10)
    
    #Generate distribution plot with genecount sums for all samples
    distribution_plot(match_id_table,rep_patient_uniqueID, mycols, cancer, "raw data")
    ggsave(paste0(cancer,"_raw_distribution_",sample_type,".png"),width = 12, height = 10)
    
    #run PCA analysis and generate plots
    MDSPlot(my.data = genecount,
            match_table = match_id_table,
            my.labels = rep_patient_uniqueID,
            my.cols = mycols,
            cancertype = as.character(cancer),
            datatype = "raw data"
    )
    ggsave(paste0(cancer,"_raw_PCA_plot.png"),width = 12, height = 8)
    
    
  }else{   
    print(paste("Replicates for", sample_type, "has not been detected"))
    discard_samples <- character(0)
  }
  
  return (discard_samples)
} 


This direcotry (Technical_replicates) perform replicate analysis compatible for running for all cancer types, one at a time.
The outputs have been ordered in folders so all the .csv files can be found for all cancer types in the replicate_info folder
and the .png files can be found for cancer types having replicates in the following folders respectively;
Barplot, Distribution_plot and PCA_plot.

* REQUIREMENTS
- Packages (version used to produce this data):
        - TCGAbiolinks (2.12.5)
        - SummarizedExperiment (1.14.1)
        - TCGAutils (1.4.0)
        - ggplot2 (3.3.5)
        - ggrepel (0.8.1)
        - limma (3.40.6)
        

* WHAT HAS BEEN DONE
Step 1. Read in RangedSummarizedExperiment object created by
Collecting female barcodes from the 'TCGAbiolinks' function 'GDCquery_clinic'
and get TCGA mRNA expression in primary tumor and normal samples from
GDC for those using TCGAbiolinks function 'GDCquery' with options:
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "HTSeq - Counts",
        sample.type = c("Primary solid Tumor",  "Solid Tissue Normal"),
        legacy = FALSE
        barcode = female_barcodes
and then aggregating the results with 'GDCprepare' using
directory = “../../GDCdata”.

Step 1 was done in the original script, but because of GDC update, data are instead
loaded, which have already been aggregated by GDCprepare from directory = "../../../data"

Step 2: Get information about replicates in data by the function check_for_replicates(),
which generate informative csv files

Step 3: Generate txt logfile describing the number of unique and total replicates

Step 4: Analyze the replicates by producing plots with the function replicate_analysis().

* OUTPUT
1. [cancer]_duplicate_info_NT.csv, [cancer]_duplicate_info_TP.csv and [cancer]_logfile.txt:
These files come from the second step where replicate information is extracted.
The .csv files are devided into NT and TP containing TCGA_barcode, TCGA_submitterID
and Match_id (refering to the replicate). The .txt file is a logfile describing the
number of unique and total number of replicates in NT and TP, as well if none 
where found.

2. [cancer]_raw_barplot.png, [cancer]_raw_distribution.png and [cancer]_raw_PCA_plot.png:
These files come from the fourth step with analyzing the replicates. 
The barplot shows gene count sums for the replicates color-coded for each patient.
The distribution plot shows the gene count sums as a scatterplot for all samples
highlighting the replicates.
The PCA plot shows all samples with the replicates highlighted.


* BUILDING THE DATABASE
There are two available R scripts: Replicate_investigation.r and
Replicate_functions.r. The Replicate_investigation is the script to run,
and it loads the functions found in Replicate_functions.r.

Run TCGA_expression.r on the terminal like:
“Rscript Replicate_functions.r [cancer]”


* INTRODUCTION *
This directory (Casestudy) performs replicate removal and pre-processing of TCGA mRNA expression data for Primary Tumor (TP) and Normal Tissue (NT) samples for Utherine Corpus Endometrial Carcinoma (UCEC). The directory contains several subfolders with individual README files explaining their purpose.

* REQUIREMENTS
- Packages (version used to produce this data):
        - TCGAbiolinks (2.12.5)
        - SummarizedExperiment (1.14.1)
        - TCGAutils (1.4.0)
        - ggplot2 (3.3.5)
        - ggrepel (0.8.1)
        - limma (3.40.6)
        

* WHAT HAS BEEN DONE *
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

Step 2: Remove samples with gene count sum < 20 million reads

Step 3: Get information about replicates in data

Step 4: Analyze the replicates by producing plots and keep only the replicate sample
with the highest gene count sum, the other samples are discarded from the data.

Step 5. Preprocessing data in three steps:
        1. Data pre-processing, sorting out samples with a spearman correlation
        coefficient less than 0.6, using function 'TCGAanalyze_Preprocessing'
        with options:
                - cor.cut = 0.6
                
        2. Data withinLane GC-count normalization using the function
        'TCGAanalyze_Normalization' with new geneInfoHT table with options:
                - geneInfo = geneInfoHT,
                - method = "gcContent"
                
        3. Data 0.25 quantile filtering with function 'TCGAanalyze_Filtering'
        with options:
                - method = "quantile",
                - qnt.cut =  0.25
        Files form all steps are saved twice, one file with ensembl IDs as row names
        and the other with HUGO names as row names.

* OUTPUT *
1. UCEC_duplicate_info_NT.csv, UCEC_duplicate_info_TP.csv and UCEC_logfile.txt:
These files come from the third step where replicate information is extracted.
The .csv files are devided into NT and TP containing TCGA_barcode, TCGA_submitterID
and Match_id (refering to the replicate). The .txt file is a logfile describing the
number of unique and total number of replicates in NT and TP, as well if none 
where found.

2. UCEC_raw_barplot.png, UCEC_raw_distribution.png and UCEC_raw_PCA_plot.png: These files
come from the fourth step with analyzing the replicates. 
The barplot shows gene count sums for the replicates color-coded for each patient.
The distribution plot shows the gene count sums as a scatterplot for all samples
highlighting the replicates.
The PCA plot shows all samples with the replicates highlighted.


3. UCEC_dataPrep.rda and UCEC_dataPrep_HUGO.rda: These files come from
the first pre-processing step and contain mRNA expression values, with TCGA
barcodes as column names and gene ID as row names: ensembl IDs in
UCEC_dataPrep.rda and HUGO names in UCEC_dataPrep_HUGO.rda
(see Step 5.1 in ‘WHAT HAS BEEN DONE’)

3. UCEC_dataNorm.rda and UCEC_dataNorm_HUGO.rda: These files come
from the GC-count normalization step and contain mRNA expression values, with
TCGA barcodes as column names and gene ID as row names: ensembl IDs
in UCEC_dataNorm.rda and HUGO names in UCEC_dataNorm_HUGO.rda
(see Step 2.2 in ‘WHAT HAS BEEN DONE’)

4. UCEC_dataFilt.rda and UCEC_dataFilt_HUGO.rda: These files come
from the quantile filtering step and contain mRNA expression values, with
TCGA barcodes as column names and gene ID as row names: ensembl IDs
in UCEC_dataFilt.rda and HUGO names in UCEC_dataNorm_HUGO.rda
(see Step 2.3 in ‘WHAT HAS BEEN DONE’)

5. dataPrep_UCEC_array_array.png: Standard output plot form the first
pre-processing step. It visualizes intensity correlation, and samples by
samples correlation.

* BUILDING THE DATABASE *
There are three available R scripts: TCGA_expression.r,
TCGA_expression_functions.r and TCGA_replicate_function.r. The
TCGA_expression.r is the script to run, and it loads the functions found in
TCGA_expression_functions.r and TCGA_replicate_function.r.

Run TCGA_expression.r on the terminal like:
“Rscript TCGA_expression.r [cancer]”

This exact line was used to produce the results found in this folder:
“Rscript TCGA_expression.r UCEC”

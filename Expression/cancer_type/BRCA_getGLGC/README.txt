This directory (BRCA_getGLGC) performs pre-processing of TCGA mRNA expression data in Primary Tumor and Normal Tissue samples
for Breast Invasive Carcinoma (BRCA). The normalization step are performed by both the new generated table and the orgiginal
geneInfoHT table to investigate the differences.

* REQUIREMENTS
- Packages (version used to produce this data):
	- TCGAbiolinks (2.12.5)
	- SummarizedExperiment (1.14.1)

* WHAT HAS BEEN DONE
Step 1. Get TCGA-BRCA mRNA expression data from GDC:
Collect female barcodes from the 'TCGAbiolinks' function 'GDCquery_clinic'
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

Step 1 was done in the original script, but because of GDC update data are instead
loaded, which have already been aggregated by GDCprepare from directory = "../../../data"

Step 2. Preprocessing data in three steps:
	1. Data pre-processing, sorting out samples with a spearman correlation 
	coefficient less than 0.6, using function 'TCGAanalyze_Preprocessing' 
	with options:
		- cor.cut = 0.6
		
	2.1 Data withinLane GC-count normalization using the function 
	‘'TCGAanalyze_Normalization' with original options:
		- geneInfo = geneInfoHT,
		- method = "gcContent"
	2.2 Data withinLane GC-count normalization using the function 
	‘'TCGAanalyze_Normalization' with new table geneinfo_GLGC:
		- geneInfo = geneinfo_GLGC,
		- method = "gcContent"
		
	3. Data 0.25 quantile filtering with function 'TCGAanalyze_Filtering' 
	with options:
		- method = "quantile",
		- qnt.cut =  0.25
	Files form all steps are saved twice, one file with ensembl IDs as row names
	and the other with HUGO names as row names.  

* OUTPUT
1. data: This is a directory containing the aggregated data downloaded from GDC
in a SummarizedExperiment format.

2. BRCA_dataPrep.rda and BRCA_dataPrep_HUGO.rda: These files come from
the first pre-processing step and contain mRNA expression values, with TCGA 
barcodes as column names and gene ID as row names: ensembl IDs in 
BRCA_dataPrep.rda and HUGO names in BRCA_dataPrep_HUGO.rda 
(see Step 2.1 in ‘WHAT HAS BEEN DONE’)

3. BRCA_dataNorm.rda and BRCA_dataNorm_HUGO.rda: These files come
from the GC-count normalization step and contain mRNA expression values, with 
TCGA barcodes as column names and gene ID as row names: ensembl IDs 
in BRCA_dataNorm.rda and HUGO names in BRCA_dataNorm_HUGO.rda 
(see Step 2.2 in ‘WHAT HAS BEEN DONE’)

4. BRCA_dataFilt.rda and BRCA_dataFilt_HUGO.rda: These files come
from the quantile filtering step and contain mRNA expression values, with 
TCGA barcodes as column names and gene ID as row names: ensembl IDs 
in BRCA_dataFilt.rda and HUGO names in BRCA_dataNorm_HUGO.rda 
(see Step 2.3 in ‘WHAT HAS BEEN DONE’)

5. dataPrep_BRCA_array_array.png: Standard output plot form the first
pre-processing step. It visualizes intensity correlation, and samples by 
samples correlation.

* BUILDING THE DATABASE
There are two available R scripts: TCGA_expression.r and 
TCGA_expression_functions.r. The TCGA_expression.r is the script to
run, and it loads the functions found in TCGA_expression_functions.r

Run TCGA_expression.r on the terminal like:
“Rscript TCGA_expression.r [cancer]”

This exact line was used to produce the results found in this folder:
“Rscript TCGA_expression.r BRCA”


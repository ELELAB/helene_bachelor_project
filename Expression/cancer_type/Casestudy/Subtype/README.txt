* INTRODUCTION *
This directory (Subtype) do differential expression analysis (DEA) by the revised TCGAanalyze function on preprocessed UCEC subtype specific data.
It compare the differential expresed genes (DEGs) by Enrichment analysis,
UpSet plots and finds the most Up- and downregulated DEGs for each subtype.

* REQUIREMENTS *
Packages (version used to produce this data):
Get_subtypes.R
  - TCGAbiolinks (2.12.5)
  - TCGAutils (1.4.0)

DEA_subtype.R
	- TCGAbiolinks (2.12.5)
	- limma (3.40.6)
	- edgeR (3.26.7)

Comparison_subtypes.R
	- grid (3.6.2)
	- gridBase (0.4-7)
	- enrichR (3.0)
	- ComplexHeatmap (2.0.0)

* WHAT HAS BEEN DONE *
Get_subtypes.R
  Step 1. UCEC pre-processed data are loaded
  
  Step 2. Data are extracted for all four subtypes (CN_LOW, CN_HIGH, MSI and POLE) 
  from the pre-processed UCEC gene count matrix and saved.
  
  Step 3. A dataframe with the number of samples in each subtypes are created
  
DEA_subtype.R
	Step 1. UCEC and subtype specific preprocessed data are loaded.

	Step 2. DEA analysis for tumor and normal sample by limma and edgeR for each of
	the subtypes are performed
	
	Step 3. Only the gene IDs which are found as DE by both limma and edgeR are kept
	in the DEA output files and saved as rda files

Comparison_subtypes.R
	Step 1. DEA output files with limma values for each subtype are loaded

	Step 2. The genes are sorted in up- and downregulating genes.
	
	Step 3. Enrichment analysis are performed for each subtype divided in up
	and downregulated DEGs against the database; 	GO_Biological_Process_2021 and
	results are plotted
	
	Step 4. Upset plots are created, visualizing intersections of DEGs found across
	subtypes and as well intersections found of up and downregulated DEGs across
	subtypes.
	
	Step 5. Genes found as both up and downregulated in the same subtype are found.
	
	Step 6. The most Up- and downregulated DEGs in each subtype are found and
	investigated in how many subtypes the DEG are found.
	

* OUTPUT *
Get_subtypes.R
  1. CN_LOW.rda, CN_HIGH.rda, MSI.rda and POLE.rda, containing subtype specific
  gene expression matrix.

DEA_subtype.R
	1. DEG_CN_LOW_Limmaboth.rda and DEG_CN_LOW_edgeRboth.rda; a dataframe with
	DEGs found for CN low and belonging values obtained from Limma or edgeR
	
	2. DEG_CN_HIGH_Limmaboth.rda and DEG_CN_HIGH_edgeRboth.rda; a dataframe with
	DEGs found for CN high and belonging values obtained from Limma or edgeR
	
	3. DEG_MSI_Limmaboth.rda and DEG_MSI_edgeRboth.rda; a dataframe with
	DEGs found for MSI and belonging values obtained from Limma or edgeR
	
	4. DEG_POLE_Limmaboth.rda and DEG_POLE_edgeRboth.rda; a dataframe with
	DEGs found for POLE and belonging values obtained from Limma or edgeR

Comparison_subtypes.R
	1. Enrichment_Limma_CN_LOW_pos.png and Enrichment_Limma_CN_LOW_neg.png;
	Enrichment barplot for CN low up or downregulated DEGs
	
	2. Enrichment_Limma_CN_HIGH_pos.png and Enrichment_Limma_CN_HIGH_neg.png;
	Enrichment barplot for CN high up or downregulated DEGs
	
	3. Enrichment_Limma_MSI_pos.png and Enrichment_Limma_MSI_neg.png;
	Enrichment barplot for MSI up or downregulated DEGs
	
	4. Enrichment_Limma_POLE_pos.png and Enrichment_Limma_POLE_neg.png;
	Enrichment barplot for POLE up or downregulated DEGs
	
	5. Upset_plot.png; Upset plot with intersections across subtypes
	
	6. Upset_UpDown_plot.png; Upset plot with intersections across subtype up and
	down regulated DEGs 

* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript DEA.analysis.R
Rscript Comparison.analysis.R

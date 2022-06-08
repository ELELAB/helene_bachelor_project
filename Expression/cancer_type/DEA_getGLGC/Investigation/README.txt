* INTRODUCTION *
This directory (Investigation) creates scatterplots with the differentially expressed genes logFC,
to compare the old and new method of normalization. Furthermore investigates the potential unique gene IDs
observed in the Venndiagram, by creating bar plots visualizing the number of potential unique genes and genes
found in the opposite table of what they were normalized by.

* REQUIREMENTS *
The directory needs to be a subfolder to DEA, as it loads data
output from the DEA.

Packages (version used to produce this data):
Unique_IDs.R
	- TCGAbiolinks (2.12.5)
	- SummarizedExperiment (1.14.1)
	- ggplot2 (3.3.5)


* WHAT HAS BEEN DONE
Non_intersection_ID.R
	Step 1. Up- and downregulated DEG are loaded for both new and old method

	Step 2. Unique IDs for up- and downregulated DEGs for the two methods are found respectively
	
	Step 3. logFC are plotted for up- and downregulated genes for both methods

Unique_IDs.R
  Step 1. Loads RangedSummarizedExperiment, new and old normalization tables

	Step 2. Finds the difference of genes present in the two tables
	
	Step 3. Loads Unique IDs for the two methods
	
	Step 4. HUGO ID is converted to ensembl to compare with the tables
	
	Step 5. Investigate if unique DEGs are present in both tables
	
	Step 6. Create dataframe for each method and visualize by bar plots
	
* OUTPUT
Non_intersection_ID.R
	1. dif_upreg_oldGCn.RDS, dif_downreg_oldGCn.RDS, dif_upreg.RDS, dif_downreg.RDS: List of unique DEGs
	2. logFC_unique_oldGCn.png and logFC_unique_new.png: scatterplot showing distribution of logFC for unique DEGs (old and new), with logFC cutoff and mean logFC marked

Unique_IDs.R
1. Barplot_unique_IDs.png: Barplot showing the number of unique new genes present in the two tables
2. Barplot_unique_IDs_oldGCn.png: Barplot showing the number of unique old genes present in the two tables

* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript Non_intersection_ID.R
Rscript Unique_IDs.R

	
	

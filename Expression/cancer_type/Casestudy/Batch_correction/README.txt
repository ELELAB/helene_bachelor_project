* INTRODUCTION *
This directory Batch_correction performs batch correction by TCGAbatch_corection on preprocessed UCEC data and visualize it in PCA plots. 

* REQUIREMENTS *
The directory needs to be subfolder to the Casestudy folder and at the same level as
the Subtype folder. UCEC and subtype specific data are read in from these two folders.

Packages (version used to produce this data):
	- TCGAbiolinks (2.12.5)
	- limma (3.40.6)
	- ggplot2 (3.3.6)

* WHAT HAS BEEN DONE
Step 1. UCEC preprocessed and subtype specific gene count matrices is loaded.

Step 2. log transformation with voom function are performed

Step 3. PCA plot of non batch corrected data and plate batch corrected data are plotted and saved. While TSS, Portion, Center and Year cannot bee used as batch factors
	
Step 4. Examine the batch factors clinical data by the summary function

* OUTPUT
1. UCEC_non_batch_corrected.png
2. UCEC_plate_batch_corrected.png

* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript Batch_correction.R
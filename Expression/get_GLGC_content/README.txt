* INTRODUCTION *
This directory (get_GLGC_content) generates a gene annotation table for GRCh38 (hg38), containing gene length and GC-content.
This directory contains as well a copy of the orginal geneInfoHT table: oldGC_geneinfoHT.txt.

* REQUIREMENTS *
Packages (version used to produce this data):
	- biomaRt (2.40.4)
	- EDASeq (2.18.0)

* WHAT HAS BEEN DONE
Step 1. ensemble_gene_id are found for hg38 by getBM function

Step 2. gene length and GC-content are found by getGeneLengthAndGCContent function for all ensemble gene id

Step 3. dataframe is created containing genelength and GC-content as columns, and ensemble gene ID as rownames. This table is saved as a rda file.

* OUTPUT
1. getGLGC_download.rda: This are the downloaded data from biomart
2. geneInfoHT.rda: This is the final file with the column names:
geneLength and gcContent, and ensemble gene ID as rownames.
3. GLGC_table.txt: This is also the final file, generated for comparison analysis with the original geneInfoHT table


* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript create_GLGC_table.r

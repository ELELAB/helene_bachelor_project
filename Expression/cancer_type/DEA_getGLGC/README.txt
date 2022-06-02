* INTRODUCTION *
This directory (DEA_getGLGC) do differential expression analysis (DEA) by the revised TCGAanalyze
function on preprocessed breast cancer data, which is normalized by a new and old
GC content normalization. It compare the differential expresed genes (DEG's) by Venndiagrams and
investigate the top 20 most up- and downregulated genes with correlations plots. It performs an enrichment analysis of DEGs from the new and old normalization method.

* REQUIREMENTS *
The directory needs to be at the same level as the BRCA folder, where the data are
generated.

Packages (version used to produce this data):
DEA.analysis.R
	- TCGAbiolinks (2.12.5)
	- limma (3.40.6)

Comparison.analysis.R
	- grid
	- futile.logger (1.4.3)
	- VennDiagram (1.7.1)

* WHAT HAS BEEN DONE
DEA.analysis.R
	Step 1. Preprocessed data by the two methods is loaded from BRCA, each in both
	ENSEMBL and HUGO ID. 

	Step 2. DEA analysis for tumor and normal sample for each of the four cases,
	are performed with fdr.cut = 0.05, logFC.cut = 0.5, pipeline = "limma",
	voom = TRUE, batch.factors ="TSS" and output are saved as rda files.

Comparison.analysis.R
	Step 1. DEG's are loaded (four files in total).

	Step 2. The genes are sorted in up-(logFC > 0) and downregulating genes (logFC < 0).
	
	Step 3. Venndiagram comparing the two methods are created and saved.
	
	Step 4. The 20 most up- and downregulated genes are tested if they are comparable
	between the two methods the data was preprocessed. An output file is saved.
	
	Step 5. Correlations plots are generated for the 20 most up- and downregulated genes.
	The plots are saved.
	
	Step 6. Enrichment analysis is performed with enrichr and the enriched biological processes plots are saved as png.


* OUTPUT
DEA.analysis.R
	1. DEG_ensembl.rda: DEGs by new GC normalizationin in ENSEMBL ID
	2. DEG_ensembl_oldGCn.rda: DEGs by old GC normalization in ENSEMBL ID
	3. DEG_HUGO.rda: DEGs by new GC normalizationin in HUGO ID
	4. DEG_HUGO_oldGCn.rda: DEGs by old GC normalization in HUGO ID

Comparison.analysis.R
	1. venn_ensembl.jpg: Venndiagram comparing the two methods in ENSEMBL ID
	2. venn_hugo.jpg: Venndiagram comparing the two methods in HUGO ID
	3. ID_comparison.txt: results from comparing the 20 most up- and downreguated
	between the two methods
	4. correlationplot_ENSEMBL.jpg: correlation plot from 20 most up- and downreguated
	between the two methods with ENSEMBL ID 
	5. correlationplot_HUGO.jpg: correlation plot from 20 most up- and downreguated
	between the two methods with HUGO ID
	6. Enrichment_hugo_down.png, Enrichment_hugo_up.png, Enrichment_hugo_oldGCn_down.png,
	Enrichment_hugo_oldGCn_up.png: Enrichmentplots showing enriched terms, genecount and p-value
	for the new DEGs up and new DEGs down, old DEGs up and old DEGs down.

* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript DEA.analysis.R
Rscript Comparison.analysis.R

* INTRODUCTION *
This directory (DEA) do differential expression analysis (DEA) by the revised TCGAanalyze
function on preprocessed UCEC data, which is normalized by the new GC content normalization.
It compare the differential expresed genes (DEG's) from limma and edgeR by Venndiagrams and
investigate the top 200 most up- and downregulated genes with correlations plots.
It performs an enrichment analysis of DEGs found by both Limma and edgeR.

* REQUIREMENTS *
Packages (version used to produce this data):
DEA.analysis.R
	- TCGAbiolinks (2.12.5)
	- limma (3.40.6)
	- edgeR (3.26.7)

Comparison.analysis.R
	- grid (3.6.2)
	- futile.logger (1.4.3)
	- VennDiagram (1.7.1)

* WHAT HAS BEEN DONE
DEA.analysis.R
	Step 1. Preprocessed data (dataFilt) are loaded with	HUGO ID. 
	
	Step 2. Rows in dataFilt without HUGO gene ID are removed and the file is saved

	Step 2. DEA analysis are performed with TCGAanalyze_DEA_revised for tumor and
	normal sample with Limma: fdr.cut = 0.1, logFC.cut = 1, pipeline = "limma",
  voom = TRUE and batch.factors ="TSS" or with edgeR: fdr.cut = 0.1, logFC.cut = 1,
  pipeline = "edgeR", method = "glmLRT", voom = FALSE and batch.factors ="TSS"
	and output is saved as rda files

Comparison.analysis.R
	Step 1. DEG's are loaded (two files in total).

	Step 2. The genes are sorted in up- (logFC > 0) and downregulating (logFC < 0) genes and saved.
	
	Step 3. Venndiagram comparing the two methods are created and saved.
	
	Step 4. Correlations plots are generated for the 200 most up- and downregulated genes,
	comparing logFC for the two methods. The plots are saved.
	
	Step 5. Enrichment analysis is performed with enrichr for the DEGs found by both limma
	and edgeR against the following databases: GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
	"GO_Biological_Process_2021","KEGG_2021_Human", "Reactome_2016","MSigDB_Hallmark_2020.
	The plots are saved as png.


* OUTPUT
DEA.analysis.R
  1. UCEC_dataFilt_HUGO_wID.rda: Gene ekspression matrix without blank gene IDs
	2. DEG_Limma.rda: DEGs found by Limma
	3. DEG_edgeR.rda: DEGs found by edgeR

Comparison.analysis.R
  1. Downreg_edgeR.RDS, Downreg_Limma.RDS, Upreg_edgeR.RDS, Upreg_Limma.RDS: 
     DEA output divided in up- and downregulated genes from limma and EdgeR
	2. venndiagram.png: Venndiagram comparing the two methods 
	3. correlationplot.png: correlation plot from 200 most up- and downreguated
	   between the two methods
	4. enriched_neg.rda, enriched_pos.rda: List of enrichment result tables for all databases
	5. Enrichment_BP_neg.png, Enrichment_BP_pos.png, Enrichment_CC_neg.png, Enrichment_CC_pos.png,
     Enrichment_Hallmark_neg.png, Enrichment_Hallmark_pos.png, Enrichment_KEGG_neg.png,
     Enrichment_KEGG_pos.png, Enrichment_MF_neg.png, Enrichment_MF_pos.png, Enrichment_Re_neg.png,
     Enrichment_Re_pos.png: Enrichmentplots show enriched terms, genecount and p-value.

* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript DEA.analysis.R UCEC
Rscript Comparison.analysis.R

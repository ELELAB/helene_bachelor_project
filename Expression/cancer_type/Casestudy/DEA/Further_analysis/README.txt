* INTRODUCTION *
This directory (Further_analysis) do visualizations (Volcano and Ridgeline plots) of the differential expressed genes (DEGs) found in the folder DEA.

* REQUIREMENTS *
The directory needs to be a subfolder to the DEA folder, where the data are
generated.

Packages (version used to produce this data):
	- TCGAbiolinks (2.12.5)
	- ggplot2 (3.3.6)
	- ggridges (0.5.3)

* WHAT HAS BEEN DONE
Step 1. DEGs found by limma and edgeR are loaded.

Step 2. DEGs found by both methods are kept, thus creating two matrices (DEGs_Limmaboth and DEGs_edgeRboth) containing the same HUGO IDs but having different values found by the respective method.

Step 3. Volcano plots are created for DEGs_Limmaboth and DEGs_edgeRboth and saved.

Step 4. Enrichment data for up and downregulated DEGs are loaded.

Step 5. The five most enriched processes for biological process, KEGG and Reactome
are found for up and downregulated DEGs respecitively.

Step 6. Ridgeline plots for the five most enriched processes, using logFC as x-axis and significance (-log10(adj.P.val)) as x-axis are created and saved.

* OUTPUT
1. volcano_Limmaboth.png
2. volcano_edgeRboth.png
3. Ridgeline_logFC.png
4. Ridgeline_adj.P.val.png

* RUNNING THE SCRIPT *
The script are run in the following way from the terminal:
Rscript Visualization_DEGs.r

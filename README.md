# Single-cell_analysis
## Introduction
In this project a public dataset (GSE183206) with information from WT and rd10 mice from Gene Expression Omnibus (GEO) repository was used.
UMAP clustering and ligand-receptor interaction analyses were carried out to characterize altered signaling dynamics and transcriptional shifts specific to different cell types. The results showed changes of intercellular communications in the diseased retina. Specifically, two critical pathological axes were identified: a metabolic support hub via the PSAP pathway, targeting the Müller glia, and a neuroinflammatory hub through the APP pathway, connecting the Microglia.

## Contents
* README.md

File with all the information about the documents and scripts contained in the repository, including the instructions on how to work with the script files.
* project_sc
    - figures
      
    Folder with UMAPs comparing the disease and healthy genotypes, and plots of pathway networks.

    - project_scripts

    **data_analysis.R**: R script with the dimensionality reduction and UMAP.
  
    **pathway_network.R**: R script for the pathway networks.

## Use
* **data_analysis.R**
    - Required packages: Seurat, SeuratDisk and ggplot2.
    - Required input: .h5 file.
    - Output: UMAP comparing genotypes side by side, UMAP comparing genotypes by overlapping them, UMAP with cell types side by side, UMAP with a specific gene highlighted, and UMAPs per cell type with genotypes overlapped.
      
Once the script file has been downloaded, the following changes could be made to the script file.
 
    - line 6: to change the name of the file. Eg. counts <- Read10X_h5("C:/Users/data.h5')
    - from line 59 to 126: to add or remove genes/cell types.
    - line 186: to choose which gene to highlight in the UMAP.

* **pathway_network.R**
    - Required packages: Seurat and CellChat.
    - Required input: .h5 file.
    - Output: chart with routes that "begin" in RD10 or are lost in WT, graph with all the detected interactions between cell type x and cell type y, figure with the pathway network of interest (rd10 and WT).
  
Once the script files has been downloaded, the following changes could be made to the script file.

    - line 9: to change the name of the file. Eg. counts <- Read10X_h5("C:/Users/data.h5')
    - from line 28 to 72: to add or remove genes/cell types.
    - line 122: to change the name of the cell clusters.
    - line 175 and 176: to change the cell types to see the interactions between. Eg. p1 <- netVisual_bubble(cellchat_WT_full, sources.use = "Muller glia", targets.use = "Microglia", title.name = "WT")
    - line 183: to filter only the populations to compare. Eg. seurat_filtrado <- subset(seurat_obj, idents = c("Muller glia", "Microglia"))
    - line 209: to see all the detected interactions between cell type x and cell type y. Eg. test_interacciones <- subsetCommunication(cellchat_RD10_full,sources.use = "Muller glia",targets.use = "Microglia")
    - line 243: to choose pathway. Eg. pathway_interes <- "APP"



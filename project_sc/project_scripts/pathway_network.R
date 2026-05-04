# ===
# Part 2 signaling pathwat network
# ===

library(Seurat)
library(CellChat)

# 1. Read the matrix from the .h5 file
counts <- Read10X_h5("C:/Users/.../GSE183206_aggr_filtered_counts_matrix.h5")

# 2. Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = "Retina_Mouse")

# 3. Extracting metadata from barcodes
seurat_obj$sample_id <- gsub(".*[-_]", "", Cells(seurat_obj))
seurat_obj$genotype <- ifelse(seurat_obj$sample_id %in% c("1", "2"), "rd10", "WildType")

# 4. Processing (Normalization and UMAP)
seurat_obj <- NormalizeData(seurat_obj) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5)

# 5. Define retina markers
markers_to_check <- list(
  "Rod photoreceptor cells" = c("Rho","Nrl","Nr2e3","Gnat1","Gnb1","Gngt1","Pde6a","Pde6b","Pde6g","Cnga1","Cngb1","Guca1a","Guca1b","Rom1","Prph2","Slc24a1","Sag","Crx","Rp1","Rp1l1","Rdh12","Abca4","Prom1","Reep6","Unc119","Tulp1","Aipl1","Elovl4","Kcnj14","Rgs9","Rgs9bp"),
  
  "Cone photoreceptor cells" = c("Opn1sw","Opn1mw","Arr3","Gnat2","Gngt2","Pde6c","Pde6h","Cnga3","Cngb3","Gnb3","Guca1a","Guca1b","Slc24a2","Slc24a4","Rxrg","Thrb","Ccdc136","Rbp3","Elovl2","Elovl5","Elovl6","Cplx3","Scg3","Olfm1"),
  
  "Rod bipolar cells" = c("Pcp2","Trpm1","Pcp4","Car8","Prkca","Vsx2","Grm6","Gpr179","Cabp5","Isl1","Ndnf","Gng13","Cacna2d3","Vstm2b","Sebox","Trnp1","Qpct","Zbtb20","Tmem215","Rgs11"),
  
  "Cone bipolar cells" = c("Vsx2","Otx2","Isl1","Cabp5","Car8","Scgn","Gria2","Grik1","Grm7","Grm8","Trpc1","Slc17a7","Nfib","Nhlh2","Bhlhe23","Satb2","Prdm13","Lhx4","Gnb3","Rasgrp1"),
  
  "Microglia" = c("Apoe","C1qa","C1qb","C1qc","Csf1r","P2ry12","Tmem119","Cx3cr1","Trem2","Tyrobp","Aif1","Hexb","Lyz2","Fcer1g","Cd68","Cd74","Cd52","Laptm5"),
  
  "Muller glia" = c("Glul","Rlbp1","Slc1a3","Aqp4","Vim","Gfap","Aldh1a1","Sox9","Sox2","Hes1","Hes5","Pax6","Nfia","Nfix","Id1","Id3","Cd44","Clu","Apoe","Sparc"),
  
  "GABA_Amacrine" = c("Gad1","Gad2","Slc6a1","Slc6a11","Gabra1","Gabra2","Gabrb2","Gabrb3","Gad1os","Gad2os","Calb1","Calb2","Pvalb","Reln","Chat"),
  
  "Horizontal cells" = c("Lhx1","Onecut1","Onecut2","Onecut3","Prox1","Tfap2b","Pax6","Barhl2",
                         "Syt2","Syt5","Snap25","Stx1b","Sv2a","Cadps2","Cplx4","Rab3a","Rph3a",
                         "Gria2","Gria3","Gria4","Grin1","Gabra2","Gabra3","Gabbr1","Slc32a1",
                         "Kcnj2","Kcnj3","Kcnh7","Cacna1a","Cacna1b","Cacna1h","Cacna2d1","Cacna2d3",
                         "L1cam","Cntnap1","Cntn1","Sdk2","Robo1","Robo2","Slit1","Slit2","Tenm1","Ptprt","Ptpro",
                         "Calb1","Ppp1r1b","Rgs7","Rgs11","Npy","Car2","Isl1",
                         "Stmn2","Stmn3","Map1a","Map1b","Tubb3","Ank2","Utrn",
                         "Megf10","Megf11","Gpr158","Ndrg4","Adarb1","Cux2"),
  
  "Starburst amacrine cells" = c("Chat","Slc18a3","Pax6","Tfap2a","Calb2"),
  
  "Reticulocytes" = c("Alas2","Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Hbq1a","Bpgm",
                      "Fech","Slc25a37","Slc25a39","Glrx5","Blvrb","Ftl1","Ncoa4",
                      "Spta1","Slc4a1","Gypa","Mpp1","Add1","Epb41l2",
                      "Bnip3l","Ucp2","Cox4i1","Cox6c","Atp5b","Atp5a1","Ndufa4",
                      "Gpx1","Gpx4","Cat","Prdx1","Sod1",
                      "Trim10","Ube2o","Lmo2","Tfrc","Carhsp1",
                      "Gapdh","Pkm","Aldoa","Ldha","Pgam1","Eno1"),
  
  "Ganglion cells" = c("Rbpms","Thy1","Tubb3","Pou4f1","Jam2","Cartpt","Pcdh9",
                       "Etv1","Kcng4","Calb2","Pvalb","Opn4"),
  
  "Pigment epithelium" = c("Rgr","LOC100045988","Pon1","Rpe65","Best1","Cdh3","Mitf",
                           "Rdh10","Arl6ip1","Rlbp1","Tbx5","Bmp4","F3","Rrh","Man1a",
                           "Sema3c","Vldlr","Atp1b1","Ctsd","Cspg5","Cldn2","Sulf1","Slc39a12",
                           "Loxl4","Slc1a1","Slc6a13","Car12","Iqgap2","Tgfa","Spon1","Flot2"),
  
  "Endothelial cells" = c("Pecam1","Cdh5","Vwf","Flt1","Kdr")
  
)



# CLASSIFICATION BY GENES

# A. Clear markers list
lista_limpia <- lapply(markers_to_check, function(x) {
  intersect(x, rownames(seurat_obj))
})

# B. Calculate the score for each group
# This tells Seurat how closely each cell matches the lists
for(i in 1:length(lista_limpia)) {
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(lista_limpia[[i]]),
    name = paste0("Score_", names(lista_limpia)[i])
  )
}

# C. Assign the name of the group with the highest score to each cell.
# Extract the scores
score_cols <- grep("Score_", colnames(seurat_obj@meta.data), value = TRUE)
score_data <- seurat_obj@meta.data[, score_cols]
colnames(score_data) <- names(lista_limpia)

# Identify the cell type with the highest score for each cell.
seurat_obj$best_match <- colnames(score_data)[max.col(score_data, ties.method = "first")]

# Assign each cluster the name that is most frequently repeated within it.
# This removes the "noise"
cluster_ann <- aggregate(best_match ~ seurat_clusters, data = seurat_obj@meta.data, 
                         FUN = function(x) names(which.max(table(x))))

# Create vector to rename
new_cluster_names <- cluster_ann$best_match
names(new_cluster_names) <- cluster_ann$seurat_clusters

# Ensure that the active identity is the cluster numbers
Idents(seurat_obj) <- "seurat_clusters"
# Check what names the clusters have right now
print(levels(seurat_obj))

# Apply the final names
seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)
seurat_obj$cell_type_final <- Idents(seurat_obj)


# 10. Cluster names
clusters_interes <- c("Rod photoreceptor cells", "Rod bipolar cells", "Cone photoreceptor cells", "Muller glia", "Cone bipolar cells", 
                      "Microglia", "Endothelial cells", "Horizontal cells", "GABA_Amacrine", "Reticulocytes", "Starburst amacrine cells")

# 11. Function to process CellChat

preparar_cc_mouse <- function(obj, group_name) {
  data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
  
# Since genes are usually symbols in .h5 files (Rho, C1qa), guarantee uniqueness:
  rownames(data.input) <- make.unique(rownames(data.input))
  
  cc <- createCellChat(object = data.input, meta = obj@meta.data, group.by = "cell_type_final")
  cc@DB <- CellChatDB.mouse # mouse database
  
  cc <- subsetData(cc) %>% 
    identifyOverExpressedGenes() %>% 
    identifyOverExpressedInteractions() %>% 
    computeCommunProb(type = "triMean", raw.use = TRUE) %>% 
    filterCommunication(min.cells = 10) %>% 
    computeCommunProbPathway() %>% 
    aggregateNet()
  return(cc)
}


# =========================================
#OPTION 1
# =========================================
# 12. Split by genotype using the complete object
seurat_RD10_full <- subset(seurat_obj, genotype == "rd10")
seurat_WT_full <- subset(seurat_obj, genotype == "WildType")

# 13. Process CellChat with the function we made (with all the cells)
cellchat_RD10_full <- preparar_cc_mouse(seurat_RD10_full)
cellchat_WT_full <- preparar_cc_mouse(seurat_WT_full)

# 14. Merging them now will give both parties data on almost all cell types
object.list_full <- list(WT = cellchat_WT_full, rd10 = cellchat_RD10_full)
cellchat_merged_full <- mergeCellChat(object.list_full, add.names = names(object.list_full))

# 15. Chart with which routes "originate" in RD10 or are lost in WT
rankNet(cellchat_merged_full, mode = "comparison", stacked = T, do.stat = T)


# ========================================
# OPTION 2: only two cell types
# ========================================
# NOT NECESSARY. Visualize ONLY the interaction of interest.
# Even if the object has 10 cell types, here we only ask for these two:

p1 <- netVisual_bubble(cellchat_WT_full, sources.use = "Muller glia", targets.use = "Microglia", title.name = "WT")
p2 <- netVisual_bubble(cellchat_RD10_full, sources.use = "Muller glia", targets.use = "Microglia", title.name = "RD10")

p1 + p2

# 12. Define the cell types of interest

# Filter only the populations to compare (e.g., Microglia and Rods)
seurat_filtrado <- subset(seurat_obj, idents = c("Rod photoreceptor cells", "Rod bipolar cells"))

# 13. Divide by genotype
seurat_wt_filtrado <- subset(seurat_filtrado, subset = genotype == "WildType")
seurat_rd10_filtrado <- subset(seurat_filtrado, subset = genotype == "rd10")

# 14. Make Seurat forget the clusters that are no longer there
seurat_wt_filtrado$cell_type_final <- droplevels(seurat_wt_filtrado$cell_type_final)
seurat_rd10_filtrado$cell_type_final <- droplevels(seurat_rd10_filtrado$cell_type_final)

# 15. Update the Idents just in case
Idents(seurat_wt_filtrado) <- seurat_wt_filtrado$cell_type_final
Idents(seurat_rd10_filtrado) <- seurat_rd10_filtrado$cell_type_final

# 16. Use function with data only from the two target cell types
cellchat_wt_filtrado <- preparar_cc_mouse(seurat_wt_filtrado)
cellchat_rd10_filtrado <- preparar_cc_mouse(seurat_rd10_filtrado)

# 17. Merge to see what signals are "born" in rd10 degeneration
object.list_filtrado <- list(WT = cellchat_wt_filtrado, rd10 = cellchat_rd10_filtrado)
cellchat_merged_filtrado <- mergeCellChat(object.list_filtrado, add.names = names(object.list_filtrado))

# 18. Compare the strength of the signaling pathways
rankNet(cellchat_merged_filtrado, mode = "comparison", stacked = TRUE, do.stat = TRUE)

# 19. Not necessary, but provide details. All detected interactions between cell type X and cell type Y
test_interacciones <- subsetCommunication(cellchat_RD10_full, 
                                          sources.use = "Muller glia", 
                                          targets.use = "Microglia")
print(test_interacciones)
# ===========================================

# NOT NECESSARY. Comparison of number of interactions and strength (weight)
gg1 <- compareInteractions(cellchat_merged_filtrado, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged_filtrado, group = c(1,2), measure = "weight")
gg1 + gg2




# =========================================
# Option 3: VIEW PATHWAY NETWORK
# =========================================
#Remove unnamed cell clusters
# 1. Only keep the cells that have an assigned name
seurat_final <- subset(seurat_obj, idents = clusters_interes)

# 2. Split and create CellChats from scratch with this clean object
seurat_WT_clean <- subset(seurat_final, genotype == "WildType")
seurat_RD10_clean <- subset(seurat_final, genotype == "rd10")

# 3. Use the previous function
cellchat_WT <- preparar_cc_mouse(seurat_WT_clean)
cellchat_RD10 <- preparar_cc_mouse(seurat_RD10_clean)

# 4. First, combine the clean objects into a list
object.list <- list(WT = cellchat_WT, RD10 = cellchat_RD10)

# 5. Choose pathway to view its network
pathway_interes <- "APP"

# 6. IMPORTANT: Calculate the maximum weight of BOTH objects for that route.
# This ensures that the arrow thickness means the same thing on both maps.
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathway_interes)

# 7. Configure the space: 1 row, 2 columns
par(mfrow = c(1,2), xpd = TRUE)

# 8. Chart for rd10
netVisual_aggregate(cellchat_RD10, 
                    signaling = pathway_interes, 
                    layout = "circle", 
                    edge.weight.max = weight.max[1], # Shared scale
                    signaling.name = paste(pathway_interes, "- rd10"))

# 9. Chart for WT
netVisual_aggregate(cellchat_WT, 
                    signaling = pathway_interes, 
                    layout = "circle", 
                    edge.weight.max = weight.max[1], # Shared scale
                    signaling.name = paste(pathway_interes, "- WT"))

# 6*. If there is no important communication on WT:

netVisual_aggregate(cellchat_RD10, 
                    signaling = pathway_interes, 
                    layout = "circle", 
                    signaling.name = paste(pathway_interes, "- rd10"))


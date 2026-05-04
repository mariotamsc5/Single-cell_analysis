library(Seurat)
library(SeuratDisk)
library(ggplot2)

# 1. Read the matrix from the .h5 file
counts <- Read10X_h5("C:/Users/.../GSE183206_aggr_filtered_counts_matrix.h5")

# 2. Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

head(Cells(seurat_obj))

# Default tags
table(seurat_obj$orig.ident)

# 3. Extract the sample identifier from the default barcode
# Taking what is after the last separator (- o _)
seurat_obj$sample_id <- gsub(".*[-_]", "", Cells(seurat_obj))

# 4. Assign Genotypes based on the GEO description
# rd10 = samples 1 y 2 | WildType = samples 3 y 4
seurat_obj$genotype <- ifelse(seurat_obj$sample_id %in% c("1", "2"), "rd10", "WildType")

# 5. Verify that they were divided correctly
table(seurat_obj$genotype)


# 6. Identify genes with high variability (important for detecting the disease)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 7. Scale data (center gene expression)
seurat_obj <- ScaleData(seurat_obj)

# 8. Principal Component Analysis (PCA)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# 9. Build UMAP
# 20-30 dimensions because it seems to be the standard for retina in Seurat v3/v4
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# 10. Finding neighbors (required before testing resolutions)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

# ===
# RESOLUTION TEST
# ===
# Try different levels and see how many clusters emerge
#for (res in seq(0.1, 1.2, by = 0.1)) {
 # seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  #n_clusters <- length(unique(Idents(seurat_obj)))
  #cat("Resolution:", res, "| Clusters found:", n_clusters, "\n")
#}

# 11. FINAL CHOICE of resolution
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# 12. Groups with genes to identify clusters
markers_to_check <- list(
  "Rod photoreceptor cells" = c("Rho","Nrl","Nr2e3","Gnat1","Gnb1","Gngt1","Pde6a","Pde6b","Pde6g","Cnga1","Cngb1","Guca1a","Guca1b","Rom1","Prph2","Slc24a1","Sag","Crx","Rp1","Rp1l1","Rdh12","Abca4","Prom1","Reep6","Unc119","Tulp1","Aipl1","Elovl4","Kcnj14","Rgs9","Rgs9bp"),
  
  "Cone photoreceptor cells" = c("Opn1sw","Opn1mw","Arr3","Gnat2","Gngt2","Pde6c","Pde6h","Cnga3","Cngb3","Gnb3","Guca1a","Guca1b","Slc24a2","Slc24a4","Rxrg","Thrb","Ccdc136","Rbp3","Elovl2","Elovl5","Elovl6","Cplx3","Scg3","Olfm1"),
  
  "Rod bipolar cells" = c("Pcp2","Trpm1","Pcp4","Car8","Prkca","Vsx2","Grm6","Gpr179","Cabp5","Isl1","Ndnf","Gng13","Cacna2d3","Vstm2b","Sebox","Trnp1","Qpct","Zbtb20","Tmem215","Rgs11"),
  
  "Cone bipolar cells" = c("Vsx2","Otx2","Isl1","Cabp5","Car8","Scgn","Gria2","Grik1","Grm7","Grm8","Trpc1","Slc17a7","Nfib","Nhlh2","Bhlhe23","Satb2","Prdm13","Lhx4","Gnb3","Rasgrp1"),
  
  "Microglia" = c("Apoe","C1qa","C1qb","C1qc","Csf1r","P2ry12","Tmem119","Cx3cr1","Trem2","Tyrobp","Aif1","Hexb","Lyz2","Fcer1g","Cd68","Cd74","Cd52","Laptm5"),
  
  "Muller glia" = c("Glul","Rlbp1","Slc1a3","Aqp4","Vim","Gfap","Aldh1a1","Sox9","Sox2","Hes1","Hes5","Pax6","Nfia","Nfix","Id1","Id3","Cd44","Clu","Apoe","Sparc"),
  
  "GABA_Amacrine" = c("Gad1","Gad2","Slc6a1","Slc6a11","Gabra1","Gabra2","Gabrb2","Gabrb3","Gad1os","Gad2os","Calb1","Calb2","Pvalb","Reln","Chat"),
  
  "Gly_Amacrine" = c("Slc6a9","Slc6a5","Glra1","Glra2","Glra3","Glra4","Glrb","Abat","Slc32a1","Pax6","Isl1","Prox1","Sox2","Bhlhe22","Calb2","Reln","Npy"),
  
  "Reticulocytes" = c("Alas2","Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Hbq1a","Bpgm",
    "Fech","Slc25a37","Slc25a39","Glrx5","Blvrb","Ftl1","Ncoa4",
    "Spta1","Slc4a1","Gypa","Mpp1","Add1","Epb41l2",
    "Bnip3l","Ucp2","Cox4i1","Cox6c","Atp5b","Atp5a1","Ndufa4",
    "Gpx1","Gpx4","Cat","Prdx1","Sod1",
    "Trim10","Ube2o","Lmo2","Tfrc","Carhsp1",
    "Gapdh","Pkm","Aldoa","Ldha","Pgam1","Eno1"),
  
  "VE-Pericytes" = c("Rgs5","Pdgfrb","Cspg4","Notch3","Myh11","Acta2","Tagln","Des","Myl9","Lmod1",
    "Tpm2","Tpm1","Myh9","Pdlim1","Csrp1",
    "Cdh5","Pecam1","Tek","Kdr","Flt1","Eng","Esam",
    "Egfl7","Aplnr","Nrp1","Robo4","Rasip1","Ephb4","Erg","Fli1",
    "Col4a1","Col4a2","Lama4","Lamb1","Lamc1","Nid1","Hspg2","Sparc",
    "Cavin1","Cav1","Anxa2","Vcl","Flna","Tns1",
    "Foxf2","Sox17","Klf2","Heyl",
    "Nos3","Ednra","Gja4","Kcnj8","Abcc9"),
  
  "Horizontal cells" = c("Lhx1","Onecut1","Onecut2","Onecut3","Prox1","Tfap2b","Pax6","Barhl2",
    "Syt2","Syt5","Snap25","Stx1b","Sv2a","Cadps2","Cplx4","Rab3a","Rph3a",
    "Gria2","Gria3","Gria4","Grin1","Gabra2","Gabra3","Gabbr1","Slc32a1",
    "Kcnj2","Kcnj3","Kcnh7","Cacna1a","Cacna1b","Cacna1h","Cacna2d1","Cacna2d3",
    "L1cam","Cntnap1","Cntn1","Sdk2","Robo1","Robo2","Slit1","Slit2","Tenm1","Ptprt","Ptpro",
    "Calb1","Ppp1r1b","Rgs7","Rgs11","Npy","Car2","Isl1",
    "Stmn2","Stmn3","Map1a","Map1b","Tubb3","Ank2","Utrn",
    "Megf10","Megf11","Gpr158","Ndrg4","Adarb1","Cux2"),
  
  "Astrocytes" = c("Gfap","Slc1a3","Aqp4","Aldh1l1","S100b"),
  
  "Starburst amacrine cells" = c("Chat","Slc18a3","Pax6","Tfap2a","Calb2"),
  
  "Ganglion cells" = c("Thy1","Pou4f2","Rbpms","Sncg","Nefl"),
  
  "Pigment epithelial cells" = c("Rpe65","Best1","Rdh5","Ttr","Mertk","Rgr","LOC100045988","Pon1","Rpe65","Best1","Cdh3","Mitf",
                                 "Rdh10","Arl6ip1","Rlbp1","Tbx5","Bmp4","F3","Rrh","Man1a",
                                 "Sema3c","Vldlr","Atp1b1","Ctsd","Cspg5","Cldn2","Sulf1","Slc39a12",
                                 "Loxl4","Slc1a1","Slc6a13","Car12","Iqgap2","Tgfa","Spon1","Flot2"),
  
  "Immune cells" = c("Ptprc","Cd3e","Cd4","Cd79a","Lyz2"),
  
  "Endothelial cells" = c("Pecam1","Cdh5","Vwf","Flt1","Kdr"),
  
  "Pericytes" = c("Rgs5","Kcnj8","Abcc9","Pdgfrb","Des"),
  
  "Fibroblasts" = c("Col1a1","Pdgfra","Fn1","Acta2","Fbln1"),
  
  "Extracellular matrix cells" = c("Col1a1","Col10a1","Lama2","Fn1","Tgfbr1"),
  
  "Cornea epithelial cells" = c("Krt3","Krt12","Klf6"),
  
  "Progenitor cells" = c("Pax6","Sox2","Rax","Hes1","Vsx2")
)


# 13. Clear markers list
lista_limpia <- lapply(markers_to_check, function(x) {
  intersect(x, rownames(seurat_obj))
})

# 14. Calculate the score for each group
# This tells Seurat how closely each cell matches the lists
for(i in 1:length(lista_limpia)) {
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(lista_limpia[[i]]),
    name = paste0("Score_", names(lista_limpia)[i])
  )
}

# 15. Assign the name of the group with the highest score to each cell.
# Extract the scores
score_cols <- grep("Score_", colnames(seurat_obj@meta.data), value = TRUE)
score_data <- seurat_obj@meta.data[, score_cols]
colnames(score_data) <- names(lista_limpia)

# 16. Identify the cell type with the highest score for each cell.
seurat_obj$best_match <- colnames(score_data)[max.col(score_data, ties.method = "first")]

# 17. Assign each cluster the name that is most frequently repeated within it.
# This removes the "noise"
cluster_ann <- aggregate(best_match ~ seurat_clusters, data = seurat_obj@meta.data, 
                         FUN = function(x) names(which.max(table(x))))

# 18. Create vector to rename
new_cluster_names <- cluster_ann$best_match
names(new_cluster_names) <- cluster_ann$seurat_clusters

# 19. Ensure that the active identity is the cluster numbers
Idents(seurat_obj) <- "seurat_clusters"

# 20. Check what names the clusters have right now
print(levels(seurat_obj))

# 21. Apply final names
seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)
seurat_obj$cell_type_final <- Idents(seurat_obj)

# 22. Final UMAP chart
DimPlot(seurat_obj, 
        reduction = "umap", 
        split.by = "genotype", 
        group.by = "cell_type_final", 
        label = TRUE, 
        label.size = 3, 
        repel = TRUE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Retina Comparative Analysis: WT vs rd10",
       subtitle = "Annotation based on Module Scoring and Clusters (Res 0.5)")

# 23. Gene Egr1 graph or another
FeaturePlot(seurat_obj, features = "Egr1", split.by = "genotype")


# 24. Define colors
colores_genotipo <- c("WildType" = "#7F7F7F", "rd10" = "#D62728")

# 25. UMAP genotypes overlap
DimPlot(seurat_obj, 
        group.by = "genotype", 
        cols = colores_genotipo,
        pt.size = 0.1) + 
  theme_minimal() +
  ggtitle("UMAP Overlap of Genotypes: WT vs rd10")

# 26. UMAP genotypes side by side
DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by = "genotype", 
        split.by = "genotype", 
        cols = colores_genotipo) + 
  theme_minimal() +
  ggtitle("WT vs rd10")

# 27. UMAP per cell type
DimPlot(seurat_obj, 
        reduction = "umap", 
        group.by = "genotype", 
        cols = colores_genotipo, 
        label = FALSE) + 
  facet_wrap(~seurat_obj$cell_type_final)

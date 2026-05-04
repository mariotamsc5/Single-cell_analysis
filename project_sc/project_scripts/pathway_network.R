# ===
# Parte 2 signaling pathwat network
# ===

library(Seurat)
library(CellChat)

# 1. Leer la matriz desde el archivo .h5
counts <- Read10X_h5("C:/Users/.../GSE183206_aggr_filtered_counts_matrix.h5")

# 2. Crear objeto Seurat
seurat_obj <- CreateSeuratObject(counts = counts, project = "Retina_Mouse")

# 3. Extraer metadatos de los Barcodes
seurat_obj$sample_id <- gsub(".*[-_]", "", Cells(seurat_obj))
seurat_obj$genotype <- ifelse(seurat_obj$sample_id %in% c("1", "2"), "rd10", "WildType")

# 4. Procesamiento (Normalización y UMAP)
seurat_obj <- NormalizeData(seurat_obj) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30) %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.5)

# 5. Definir marcadores clásicos de retina (Mouse)
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



# CLASIFICACIÓN POR GENES

# A. Limpiar lista de marcadores
lista_limpia <- lapply(markers_to_check, function(x) {
  intersect(x, rownames(seurat_obj))
})

# B. Calcular el puntaje para cada uno de los grupos
# Esto le dice a Seurat cuanto se parece cada célula a las listas
for(i in 1:length(lista_limpia)) {
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(lista_limpia[[i]]),
    name = paste0("Score_", names(lista_limpia)[i])
  )
}

# C. Asignar el nombre del grupo con el puntaje más alto a cada célula
# Extraer los scores
score_cols <- grep("Score_", colnames(seurat_obj@meta.data), value = TRUE)
score_data <- seurat_obj@meta.data[, score_cols]
colnames(score_data) <- names(lista_limpia)

# Identificar el tipo celular con score máximo para cada célula
seurat_obj$best_match <- colnames(score_data)[max.col(score_data, ties.method = "first")]

# Asignar a CADA CLUSTER el nombre que más se repite dentro de él
# Esto limpia el "ruido"
cluster_ann <- aggregate(best_match ~ seurat_clusters, data = seurat_obj@meta.data, 
                         FUN = function(x) names(which.max(table(x))))

# Crear vector para renombrar
new_cluster_names <- cluster_ann$best_match
names(new_cluster_names) <- cluster_ann$seurat_clusters

# Aseguramos que la identidad activa sean los números de los clusters
Idents(seurat_obj) <- "seurat_clusters"
# Verificamos qué nombres tienen los clusters ahora mismo
print(levels(seurat_obj))

# Aplicar los nombres finales
seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)
seurat_obj$cell_type_final <- Idents(seurat_obj)


# 10. Nombres de clusters
clusters_interes <- c("Rod photoreceptor cells", "Rod bipolar cells", "Cone photoreceptor cells", "Muller glia", "Cone bipolar cells", 
                      "Microglia", "Endothelial cells", "Horizontal cells", "GABA_Amacrine", "Reticulocytes", "Starburst amacrine cells")

# 11. Función para procesar CellChat (Ahorra repetir código)

preparar_cc_mouse <- function(obj, group_name) {
  data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
  
  # Si los genes fueran Ensembl, aquí iría el mapIds(org.Mm.eg.db...)
  # Como en .h5 ya suelen ser Símbolos (Rho, C1qa), solo aseguro unicidad:
  rownames(data.input) <- make.unique(rownames(data.input))
  
  cc <- createCellChat(object = data.input, meta = obj@meta.data, group.by = "cell_type_final")
  cc@DB <- CellChatDB.mouse # BASE DE DATOS RATÓN
  
  cc <- subsetData(cc) %>% 
    identifyOverExpressedGenes() %>% 
    identifyOverExpressedInteractions() %>% 
    computeCommunProb(type = "triMean", raw.use = TRUE) %>% 
    filterCommunication(min.cells = 10) %>% 
    computeCommunProbPathway() %>% 
    aggregateNet()
  return(cc)
}


# <-- hasta aquí para todo igual

# =========================================
#OPCION 1
# =========================================
# 12. Dividir por genotipo usando el objeto COMPLETO
seurat_RD10_full <- subset(seurat_obj, genotype == "rd10")
seurat_WT_full <- subset(seurat_obj, genotype == "WildType")

# 13. Procesar CellChat con la función que hicimos (pero con todas las células)
cellchat_RD10_full <- preparar_cc_mouse(seurat_RD10_full)
cellchat_WT_full <- preparar_cc_mouse(seurat_WT_full)

# 14. Fusionar ahora que ambos tendrán datos de casi todos los tipos celulares
object.list_full <- list(WT = cellchat_WT_full, rd10 = cellchat_RD10_full)
cellchat_merged_full <- mergeCellChat(object.list_full, add.names = names(object.list_full))

# 15. Este gráfico te dirá qué vías "nacen" en RD10 o se pierden en WT
rankNet(cellchat_merged_full, mode = "comparison", stacked = T, do.stat = T)


# ========================================
# solo dos tipos de celulas
# ========================================
# NO NECESARIO. Visualizar SÓLO la interacción que interesa
# Aunque el objeto tenga 10 tipos de células, aquí le pedimos solo estas dos:

p1 <- netVisual_bubble(cellchat_WT_full, sources.use = "Muller glia", targets.use = "Microglia", title.name = "WT")
p2 <- netVisual_bubble(cellchat_RD10_full, sources.use = "Muller glia", targets.use = "Microglia", title.name = "RD10")

p1 + p2

# 12. Definimos los tipos celulares de interés

# Filtrar solo las poblaciones que quieres comparar (ej: Microglia y Rods)
seurat_filtrado <- subset(seurat_obj, idents = c("Rod photoreceptor cells", "Rod bipolar cells"))

# 13. Dividir por genotipo
seurat_wt_filtrado <- subset(seurat_filtrado, subset = genotype == "WildType")
seurat_rd10_filtrado <- subset(seurat_filtrado, subset = genotype == "rd10")

# 14. Forzar a Seurat a olvidar los clusters que ya no están
seurat_wt_filtrado$cell_type_final <- droplevels(seurat_wt_filtrado$cell_type_final)
seurat_rd10_filtrado$cell_type_final <- droplevels(seurat_rd10_filtrado$cell_type_final)

# 15. Actualizar los Idents por si acaso
Idents(seurat_wt_filtrado) <- seurat_wt_filtrado$cell_type_final
Idents(seurat_rd10_filtrado) <- seurat_rd10_filtrado$cell_type_final

# 16. Usar función con datos solo de los dos tipos celulares target
cellchat_wt_filtrado <- preparar_cc_mouse(seurat_wt_filtrado)
cellchat_rd10_filtrado <- preparar_cc_mouse(seurat_rd10_filtrado)

# 17. Fusión para ver qué señales "nacen" en la degeneración rd10
object.list_filtrado <- list(WT = cellchat_wt_filtrado, rd10 = cellchat_rd10_filtrado)
cellchat_merged_filtrado <- mergeCellChat(object.list_filtrado, add.names = names(object.list_filtrado))

# 18. Comparar la fuerza de las vías de señalización
rankNet(cellchat_merged_filtrado, mode = "comparison", stacked = TRUE, do.stat = TRUE)

# 19. NECESARIO PERO DA DETALLES. TODAS las interacciones detectadas entre cell type x y cell type y
test_interacciones <- subsetCommunication(cellchat_RD10_full, 
                                          sources.use = "Muller glia", 
                                          targets.use = "Microglia")
print(test_interacciones)
# ===========================================

# NO NECESARIO. Comparación de número de interacciones y fuerza (weight)
gg1 <- compareInteractions(cellchat_merged_filtrado, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged_filtrado, group = c(1,2), measure = "weight")
gg1 + gg2




# =========================================
# VISUALIZAR PATHWAY NETWORK
# =========================================
#eliminar clusters sin nombre celular
# 1. Ir al objeto Seurat original y quitar los clusters con números
# Solo nos quedamos con las células que tienen un nombre asignado
seurat_final <- subset(seurat_obj, idents = clusters_interes)

# 2. Dividir y crear los CellChat desde cero con este objeto limpio
seurat_WT_clean <- subset(seurat_final, genotype == "WildType")
seurat_RD10_clean <- subset(seurat_final, genotype == "rd10")

# 3. Utilizar la función de antes
cellchat_WT <- preparar_cc_mouse(seurat_WT_clean)
cellchat_RD10 <- preparar_cc_mouse(seurat_RD10_clean)

# 4. Primero, unir los objetos limpios en una lista
object.list <- list(WT = cellchat_WT, RD10 = cellchat_RD10)

# 5. Escoger pathway para ver su network
pathway_interes <- "APP"

# 6. IMPORTANTE: Calcular el peso máximo entre AMBOS objetos para esa vía
# Esto asegura que el grosor de las flechas signifique lo mismo en los dos mapas
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathway_interes)

# 7. Configuramos el espacio: 1 fila, 2 columnas
par(mfrow = c(1,2), xpd = TRUE)

# 8. Gráfico para RD10 usando weight.max
netVisual_aggregate(cellchat_RD10, 
                    signaling = pathway_interes, 
                    layout = "circle", 
                    edge.weight.max = weight.max[1], # Escala compartida
                    signaling.name = paste(pathway_interes, "- rd10"))

# 9. Gráfico para WT usando weight.max
netVisual_aggregate(cellchat_WT, 
                    signaling = pathway_interes, 
                    layout = "circle", 
                    edge.weight.max = weight.max[1], # Escala compartida
                    signaling.name = paste(pathway_interes, "- WT"))

# 6*. Si no hay comunicación importante en wt:

netVisual_aggregate(cellchat_RD10, 
                    signaling = pathway_interes, 
                    layout = "circle", 
                    signaling.name = paste(pathway_interes, "- rd10"))




#-----------------------Detailed analysis for each cluster-------------------------

cluster1 <- subset(combined2, idents = "Cluster1")
DefaultAssay(Cluster1) <- "SCT"


# PCA
FAPs_cluster <- RunPCA(FAPs_cluster)

# Apply Harmony to correct for batch effects 
FAPs_cluster <- RunHarmony(FAPs_cluster, group.by.vars = "orig.ident")

# Run UMAP on Harmony-reduced dimensions
FAPs_cluster <- RunUMAP(FAPs_cluster, dims = 1:20)
FAPs_cluster <- FindNeighbors(FAPs_cluster, dims = 1:20)
FAPs_cluster <- FindClusters(FAPs_cluster, resolution = 0.5)


# -------------- Visualization ----------------------------------------


FeaturePlot(Cluster1, features = c("gene1", "gene2"), split.by = 'orig.ident')
VlnPlot(Cluster1, features = c("gene1", "gene2"), split.by = 'orig.ident',split.plot = TRUE, cols = c("#8c96a5", "#c994c7"), pt.size = 0) 
DotPlot(cluster1, features = c("gene1", "gene2"), split.by = "orig.ident", cols=c("lightgrey","red3")) + RotatedAxis()

DimPlot(cluster1, reduction = "umap", label = T,label.size = 3.0, pt.size = 0.2, cols = cluster_colors, split.by = "orig.ident") + plot_annotation(title = "Cluster1")



#find DEG
FAPs_cluster <- PrepSCTFindMarkers(FAPs_cluster)

FAPs_all_markers <- FindAllMarkers(
  FAPs_cluster,
  only.pos = TRUE,         # 각 클러스터에서 높게 발현된 유전자만
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "SCT",
  recorrect_umi = FALSE)

cluster_DEG <- FindMarkers(
  FAPs_cluster,
  ident.1 = 3,
  ident.2 = 7,
  assay = "SCT",
  recorrect_umi = FALSE)







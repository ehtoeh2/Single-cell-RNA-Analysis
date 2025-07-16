#Road required library
library(openxlsx)  
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(writexl)
library(harmony)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(fgsea)
library(enrichplot)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(tidyverse)
library(scales)

#Set seed for reproducibility 
set.seed(1234)  



# Load 10X Genomics data for each sample
control.data1 <- Read10X(data.dir = "/My path to/Control_filtered_feature_bc_matrix")
control.data2 <- Read10X(data.dir = "/My path to/Control2_filtered_feature_bc_matrix")
Experiment.data1 <- Read10X(data.dir = "/My path to/EXP_filtered_feature_bc_matrix")
Experiment.data2 <- Read10X(data.dir = "/My path to/EXP2_filtered_feature_bc_matrix")


# Create Seurat objects
Control1 <- CreateSeuratObject(control.data1, project = "Control")
Control2 <- CreateSeuratObject(control.data2, project = "Control")
Experiment1 <- CreateSeuratObject(Experiment.data1, project = "Experiment")
Experiment2 <- CreateSeuratObject(Experiment.data2, project = "Experiment")



# Calculate percentage of mitochondrial genes
Control1[["percent.mt"]] <- PercentageFeatureSet(Control1, pattern = "^mt-")
Control2[["percent.mt"]] <- PercentageFeatureSet(Control2, pattern = "^mt-")
Experiment1[["percent.mt"]] <- PercentageFeatureSet(Experiment1, pattern = "^mt-")
Experiment2[["percent.mt"]] <- PercentageFeatureSet(Experiment2, pattern = "^mt-")


# Calculate percentage of ribosomal genes
Control1[["percent.ribo"]] <- PercentageFeatureSet(Control1, pattern = "^Rps|^Rpl")
Control2[["percent.ribo"]] <- PercentageFeatureSet(Control2, pattern = "^Rps|^Rpl")
Experiment1[["percent.ribo"]] <- PercentageFeatureSet(Experiment1, pattern = "^Rps|^Rpl")
Experiment2[["percent.ribo"]] <- PercentageFeatureSet(Experiment2, pattern = "^Rps|^Rpl")



# Visualize QC metrics (features, counts, mt, ribo)
VlnPlot(Control1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(Control2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(Experiment1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(Experiment2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)


# Filter out low-quality cells
Control1 <- subset(Control1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) # Add ribosomal genes if you want (optional)
Control2 <- subset(Control2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
Experiment1 <- subset(Experiment1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
Experiment2 <- subset(Experiment2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)


# Merge samples into groups
Control_merged <- merge(Control1, y = Control2, add.cell.ids = c("Control1", "Control2"))
Experiment_merged <- merge(Experiment1, y = Experiment2, add.cell.ids = c("Experiment1", "Experiment2"))


# Merge Control and Experiment groups together
combined <- merge(Control_merged, y = Experiment_merged, add.cell.ids = c("Control","Experiment"))



# Run SCTransform normalization, regressing out mitochondrial content
options(future.globals.maxSize = 1e9)  # Set max memory to 1GB for future operations
DefaultAssay(combined) <- "RNA"
combined <- SCTransform(combined, vars.to.regress = "percent.mt", verbose = FALSE)



# Run PCA
combined <- RunPCA(combined, npcs = 50)


# Visualize PCA elbow plot to determine optimal number of PCs
ElbowPlot(combined, ndims = 50)


# Run Harmony for batch correction and Find neighbors and clusters 
combined <- RunHarmony(combined, group.by.vars = "orig.ident", dims = 1:20)
combined <- FindNeighbors(combined, resolution = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)


# Run UMAP on Harmony dimensions for visualization
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)


# Set orig.ident factor levels for plotting
combined$orig.ident <- factor(combined$orig.ident, levels = c("Control", "Experiment"))


# Plot UMAP with cluster labels and sample splits
DimPlot(combined, reduction = "umap", label = T, label.size = 3.0, pt.size = 0.2) + plot_annotation(title = "Integration single cell data")
DimPlot(combined, reduction = "umap", label = T, split.by = "orig.ident", label.size = 3.0, pt.size = 0.2) + plot_annotation(title = "Integration single cell data")



# Prepare object for marker gene analysis
combined <- PrepSCTFindMarkers(combined)


# Identify marker genes for each cluster
combined.marker <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)





#Visualize marker genes (Dot Plot)
DefaultAssay(combined) <- "SCT"

marker_genes <- c("Gene1","Gene2","Gene3","Gene4","Gene5","Gene6")

# Visualize expression of marker genes across cells 
FeaturePlot(combined, features = marker_genes)
DotPlot(combined, features = marker_genes, cols=c("lightgrey","red3")) + RotatedAxis()
VlnPlot(combined, features = marker_genes)
        

# Assign new cluster names (manually)
combined2 <- combined

# Define new cluster IDs (number of names must match number of clusters!)
new.cluster.ids <- c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")
names(new.cluster.ids) <- levels(combined2)

# Rename cluster identities
combined2 <- RenameIdents(combined2, new.cluster.ids)

# Set new cluster ordering for plotting consistency
new_order <- c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5")
combined2 <- SetIdent(combined2, value = factor(Idents(combined2), levels = new_order))

# Set sample order for faceting (Control → Experiment)
combined2$orig.ident <- factor(combined2$orig.ident, levels = c("Control", "Experiment"))


#------------------------Visualize again after annotation--------------------------------------


#After annotation
DotPlot(combined2, features = marker_genes, cols=c("lightgrey","red3")) + RotatedAxis()
FeaturePlot(combined2, features = marker_genes)


cluster_colors <- c(
  "Cluster1"        = "#a6cee3",  
  "Cluster2"        = "#1f78b4",  
  "Cluster3"        = "#8dd3c7",  
  "Cluster4"        = "#4daf4a",  
  "Cluster5"        = "#bebada"  
)



DimPlot(combined2, reduction = "umap", label = T, label.size = 3.0, pt.size = 0.2, cols = cluster_colors) + plot_annotation(title = "Integration single cell data")



#--------------------Stacked vln plot---------------------------------------------------------



markergenes_up <- c("gene1","gene2","gene3","gene4","gene5") 



b <- VlnPlot(combined2, markergenes_up, stack = TRUE, sort = FALSE, flip = TRUE, 
             fill.by = "ident", cols = cluster_colors, adjust = 3) +
  theme_classic()+
  theme(legend.position = "right",
        scale_y_discrete(position = "left"),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, face = "bold"),
        strip.background = element_blank(),
        strip.placement = "left"
  )+  
  ylab("") +  
  ggtitle("Upregulated DEG")

b





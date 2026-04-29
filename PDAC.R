library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(ggplot2)
library(future)
library(clustree)

#### Cluster ####

options(future.globals.maxSize = 120 * 1024^3)
plan("multicore", workers = 1)
setwd("/mnt/comics-data/comics01/Objeto.PDAC")


#------------------------------     All    -------------------------------------

#### Subset (only P1-P24) ####

patients <- c(paste0("P0", 1:9), paste0("P", 10:24))

seurat <- subset(seurat, subset = sample %in% patients)


#### QC #### 

# % mitochondrial

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, "^MT")

# QC metrics (violin plot)

p <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_VlnPlot.pdf", plot = p , width = 10, height = 8)

# filter

seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25 )


#### Assays (RNA and TE) ####

# Identify which rows are TEs and which are genes

all_features <- rownames(seurat[["RNA"]])
te_features <- all_features[grepl("^TE-", all_features)]
gene_features <- all_features[!grepl("^TE-", all_features)]

# Extract the TE count matrix
# We pull from the 'counts' slot to get the raw data
te_counts <- GetAssayData(seurat, assay = "RNA", layer = "counts")[te_features, ]

# Create a new assay for the TEs
seurat[["TE"]] <- CreateAssayObject(counts = te_counts)

# Remove the TE rows from the original RNA assay
# The cleanest way to do this is to subset the assay to only include gene_features
seurat[["RNA"]] <- subset(seurat[["RNA"]], features = gene_features)


#...................................    RNA    .................................


# Set default Assay
DefaultAssay(seurat) <- "RNA"


#### Normalization and Scaling ####

#Normalize
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
seurat <- FindVariableFeatures(object = seurat, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
seurat <- ScaleData(object = seurat, features = VariableFeatures(seurat))

#Run PCA and check for possible separations based on specific variables (percent mt)

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat))


p <- DimPlot(seurat, reduction = "pca", group.by = "percent.mt") + NoLegend()
ggsave("PCA_percentmt.pdf", plot = p , width = 10, height = 8)


#### RUN SCT ####

seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

seurat <- RunPCA(object = seurat, features = VariableFeatures(object = seurat))


#### Run Harmony and Determine Dimensions for 90% Variance ####

seurat <- RunHarmony(seurat, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- seurat@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#### Find Neighbors and Clusters ####

#Find Neighbors 

seurat <- FindNeighbors(seurat, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

#Find Clusters

seurat <- FindClusters(seurat, verbose = TRUE, resolution = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))

p <- clustree(seurat)
ggsave("Clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

# Run UMAP
seurat <- RunUMAP(seurat, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

saveRDS(seurat,"PDAC_UMAP.rds")-------------------------------------------------
  
Idents(seurat) <- "SCT_snn_res.0.2"

#UMAP clusters

p <- DimPlot(seurat, reduction = "umap")
ggsave("UMAP_clusters.pdf", plot = p , width = 10, height = 8)

#UMAP patients 

p <- DimPlot(seurat, reduction = "umap", group.by = "sample")
ggsave("UMAP_samples.pdf", plot = p , width = 10, height = 8)


#### DotPlot ####

features_to_plot <- c("CCL5","GNLY","NKG7","CCL4","IL7R","CD3D","TFF1","KRT19","LCN2","SPINK1","TFF2","PRSS1","SPARCL1","ACKR1","CCL14","VWF","COL4A1","CDH5","SPP1","APOE","S100A9","IFI30","HLA-DRA","UBE2C","AKR1C2","UBE2S","CCNB1","S100A2","HIST1H4C","HMGB2","STMN1","MKI67","CENPF","IGKC","IGLC2","IGHG1","IGLC3","IGLC1","MS4A1","MMP9","CTSK","ACP5","CST3","RGS10","AIF1","COL3A1","COL1A2","COL1A1","SFRP2","DCN","LUM","TPSAB1","TPSB2","CPA3","MS4A2","rna_KIT","RGS5","CHGB")

p <- DotPlot(object = seurat, features = features_to_plot, cols = "RdBu", 
             group.by = "SCT_snn_res.0.2", assay='SCT') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("DotPlot_datapaper.pdf"), plot = combined_plot, width = 13, height = 8)


#### Find Markers ####

seurat <- PrepSCTFindMarkers(seurat)

seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, assay = "SCT",
                                 logfc.threshold = 1,min.pct=0.2, group.by = "SCT_snn_res.0.2")

#Top 15 markers

top15_markers <- seurat.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "Top15_Markers_per_Cluster.csv", row.names = FALSE)


#### Annotation (column) ####

clusters <- FetchData(seurat, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "Epithelial",
  "1" = "T/NK",
  "2" = "Myeloid",
  "3" = "T/NK",
  "4" = "Epithelial",
  "5" = "Erythrocyte",
  "6" = "CAFs",
  "7" = "Myeloid",
  "8" = "Epithelial",
  "9" = "Enterocyte",
  "10" = "B/Plasma",
  "11" = "Endothelial",
  "12" = "B/Plasma",
  "13" = "Prolif. epithelial",
  "14" = "Epithelial",
  "15" = "CAFs",
  "16" = "Epithelial",
  "17" = "Mast",
  "18" = "B/Plasma",
  "19" = "Endocrine",
  "20" = "CAFs",
  "21" = "Epithelial",
  "22" = "Epithelial",
  "23" = "T/NK"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

seurat[["cell type"]] <- celltypes

saveRDS(seurat,"PDAC_celltypes.rds")--------------------------------------------

#UMAP grup.by cell type

Idents(seurat) <- "SCT_snn_res.0.2"  
p <- DimPlot(seurat, reduction = "umap", group.by = "cell type")
ggsave("UMAP_name_celltypes.pdf", plot = p , width = 10, height = 8)

#VlnPlots
Idents(seurat) <- "cell type"

p <- VlnPlot(seurat, features = "MS4A1")
ggsave("VlnPlot_MS4A1.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "COL3A1")
ggsave("VlnPlot_COL3A1.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "CHGB")
ggsave("VlnPlot_CHGB.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "CDH5")
ggsave("VlnPlot_CDH5.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "APOA1")
ggsave("VlnPlot_APOA1.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "KRT19")
ggsave("VlnPlot_KRT19.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "HBB")
ggsave("VlnPlot_HBB.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "KIT")
ggsave("VlnPlot_KIT.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "CD68")
ggsave("VlnPlot_CD68.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "MKI67")
ggsave("VlnPlot_MKI67.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(seurat, features = "CD3E")
ggsave("VlnPlot_CD3E.pdf", plot = p , width = 10, height = 8)

#...................................    TEs    .................................


Idents(seurat) <- "cell type"

DefaultAssay(seurat) <- "TE"

#### Normalization and Scaling ####

#Normalize
seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
seurat <- FindVariableFeatures(object = seurat, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
seurat <- ScaleData(object = seurat, features = VariableFeatures(seurat))

saveRDS(seurat,"PDAC_TE.rds")---------------------------------------------------

#### TOP 15 TEs ####

te_matrix <- GetAssayData(seurat, assay = "TE", layer = "data")

te_means <- rowMeans(as.matrix(te_matrix))

top15_tes <- sort(te_means, decreasing = TRUE)[1:15]

print(top15_tes)

#Feature Plot

top15 <- c("TE-AluSx", "TE-AluSz", "TE-AluSx1", "TE-AluY", "TE-AluJb", "TE-AluSp", 
           "TE-AluJr", "TE-L2a", "TE-AluSq2", "TE-AluSg", "TE-MIRb", "TE-AluSc", "TE-MIR",
           "TE-AluSx3", "TE-AluSz6")



for (x in top15){
  p <- FeaturePlot(seurat, features = x, reduction = "umap")
  file_name <- paste0("FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#### FindMarkers ####

#FindMarkers (annotation)

seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, assay = "TE",
                               logfc.threshold = 1,min.pct=0.2, group.by = "cell type")


top15_markers <- seurat.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "TE_Top15_Markers_per_celltype.csv", row.names = FALSE)

# Feature Plot

te <- c("TE-(T)n", "TE-MIR", "TE-L1ME1", "TE-AluSg4", "TE-AluSq",
        "TE-AluJo")



for (x in te){
  p <- FeaturePlot(seurat, features = x, reduction = "umap")
  file_name <- paste0("TEmarkers_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}


#---------------------------     Epithelial   ----------------------------------

Idents(seurat) <- "cell type"
DefaultAssay(seurat) <- "RNA"

#### Subset ####

epit <- subset(seurat, idents = "Epithelial")


#...................................    RNA    .................................


#### Normalization and Scaling ####

#Normalize
epit <- NormalizeData(object = epit, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
epit <- FindVariableFeatures(object = epit, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
epit <- ScaleData(object = epit, features = VariableFeatures(epit))


#### Run SCT ####

epit <- SCTransform(epit, vars.to.regress = "percent.mt", verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

epit <- RunPCA(object = epit, features = VariableFeatures(object = epit))


#### Run Harmony and Determine Dimensions for 90% Variance ####

epit <- RunHarmony(epit, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- epit@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#### Find Neighbors e Clusters ####

#Find Neighbors 

epit <- FindNeighbors(epit, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

#Find Clusters

epit <- FindClusters(epit, verbose = TRUE, resolution = c(0, 0.05, 0.1, 0.2, 0.4, 0.6))

p <- clustree(epit)
ggsave("Epit_clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

epit <- RunUMAP(epit, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

saveRDS(epit,"Epit_UMAP.rds")---------------------------------------------------

Idents(epit) <- "SCT_snn_res.0.2"

#UMAP clusters

p <- DimPlot(epit, reduction = "umap")
ggsave("Epit_UMAP_clusters.pdf", plot = p , width = 10, height = 8)

#UMAP patients 

p <- DimPlot(epit, reduction = "umap", group.by = "sample")
ggsave("Epit_UMAP_samples.pdf", plot = p , width = 10, height = 8)

#UMAP Moffitt

epit[["moffitt"]] <- NA
epit$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", epit$sample), "Basal",
                       ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", epit$sample), "Classical",
                              "Intermediate"))


p <- DimPlot(epit, reduction = "umap", group.by = "SCT_snn_res.0.2", split.by = "moffitt")
ggsave("Epit_UMAP_moffitt.pdf", plot = p , width = 14, height = 8)


#### DotPlot ####

features_to_plot <- c("CTRB2", "AMY2A", "PRSS1", "CA2", "SLC4A4", "CFTR", "AQP1", "MUC6", "ONECUT2", "MUC5AC", "TFF1", "CLND18", "GATA6", "AGR2", "APOA2", "KRT19")

p <- DotPlot(object = epit, features = features_to_plot, cols = "RdBu", 
             group.by = "SCT_snn_res.0.2", assay='SCT') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("Epit_DotPlot.pdf"), plot = combined_plot, width = 10, height = 8)

#### FindMarkers ####

epit <- PrepSCTFindMarkers(epit)

epit.markers <- FindAllMarkers(epit, only.pos = TRUE, assay = "SCT",
                               logfc.threshold = 1,min.pct=0.2, group.by = "SCT_snn_res.0.2")

top15_markers <- epit.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "Epit_genes_Top15_Markers_per_Cluster.csv", row.names = FALSE)

#### KRT19 and EPCAM VlnPlot ####

p <- VlnPlot(epit, features = "KRT19")
ggsave("Epit_VlnPlot_KRT19.pdf", plot = p , width = 10, height = 8)

p <- VlnPlot(epit, features = "EPCAM")
ggsave("Epit_VlnPlot_EPCAM.pdf", plot = p , width = 10, height = 8)

#### Annotation Column ####

clusters <- FetchData(epit, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "Other epithelial",
  "1" = "Tumor (stressed)",
  "2" = "T/NK",
  "3" = "PanIN",
  "4" = "ADM",
  "5" = "Tumor (basal-like)",
  "6" = "Tumor (intermediate)",
  "7" = "Tumor (classical)",
  "8" = "Tumor (basal-like)",
  "9" = "Tumor (classical)",
  "10" = "Tumor (proliferative)",
  "11" = "Erythrocytes",
  "12" = "Acinar",
  "13" = "Myeloid",
  "14" = "Tumor (stressed)",
  "15" = "Tumor (basal-like)",
  "16" = "CAFs",
  "17" = "Ductal",
  "18" = "Myeloid",
  "19" = "Tumor (intermediate)"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

epit[["annotation"]] <- celltypes

saveRDS(epit,"Epit_celltypes.rds")--------------------------------------------

#UMAP group.by cell type
  
Idents(epit) <- "SCT_snn_res.0.2"  
p <- DimPlot(epit, reduction = "umap", group.by = "annotation")
ggsave("Epit_UMAP_group.by_celltypes.pdf", plot = p , width = 10, height = 8)

#### BarPlot ####

#subset and column

epit_type <- subset(epit, subset = SCT_snn_res.0.2 %in% c("0", "1", "3", "4", "5", "6", "7", "8", "9", "10", "12", "14", "15", "17", "19"))

clusters <- FetchData(epit_type, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "Other epithelial",
  "1" = "Tumor (stressed)",
  "3" = "PanIN",
  "4" = "ADM",
  "5" = "Tumor (basal-like)",
  "6" = "Tumor (intermediate)",
  "7" = "Tumor (classical)",
  "8" = "Tumor (basal-like)",
  "9" = "Tumor (classical)",
  "10" = "Tumor (proliferative)",
  "12" = "Acinar",
  "14" = "Tumor (stressed)",
  "15" = "Tumor (basal-like)",
  "17" = "Ductal",
  "19" = "Tumor (intermediate)"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

epit_type[["Epithelial types"]] <- celltypes

#barplot

plot_data <- epit_type@meta.data %>%
  group_by(orig.ident, `Epithelial types`) %>% 
  summarise(cnt = n()) %>%
  group_by(orig.ident) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  ungroup()

p <- ggplot(plot_data, aes(x = orig.ident, y = freq, fill = `Epithelial types`)) +
  geom_bar(stat = "identity", width = 0.8) + 
  theme_bw() + 
  labs(
    y = "Proportion of cells in cluster", 
    x = NULL, 
    fill = NULL
  ) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "bottom", 
    panel.grid = element_blank(), 
    axis.line = element_line(colour = "black")
  )

ggsave("Epit_barplot.pdf", plot = p , width = 10, height = 8)#barplot
#...................................    TEs    .................................


DefaultAssay(epit) <- "TE"


#### Normalization and Scaling ####

#Normalize
epit <- NormalizeData(object = epit, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
epit <- FindVariableFeatures(object = epit, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
epit <- ScaleData(object = epit, features = VariableFeatures(epit))

saveRDS(epit,"Epit_TE.rds")-----------------------------------------------------

  
#### TOP 15 TEs ####

te_matrix <- GetAssayData(epit, assay = "TE", layer = "data")

te_means <- rowMeans(as.matrix(te_matrix))

top15_tes <- sort(te_means, decreasing = TRUE)[1:15]

print(top15_tes)


#Feature Plot

top15 <- c("TE-AluSx", "TE-AluSz", "TE-AluSx1", "TE-AluY", "TE-AluJb", "TE-AluSp", 
           "TE-AluJr", "TE-L2a", "TE-AluSq2", "TE-AluSg", "TE-MIRb", "TE-AluSc", "TE-MIR",
           "TE-AluSx3", "TE-AluSz6")



for (x in top15){
  p <- FeaturePlot(epit, features = x, reduction = "umap")
  file_name <- paste0("Epit_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}


#### FindMarkers ####

#FindMarkers (annotation)

epit.markers <- FindAllMarkers(epit, only.pos = TRUE, assay = "TE",
                               logfc.threshold = 1,min.pct=0.2, group.by = "annotation")


top15_markers <- epit.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "Epit_TE_Top15_Markers_per_celltype.csv", row.names = FALSE)

# Feature Plot

te <- c("TE-Harlequin-int", "TE-LTR48B", "TE-HERVH-int", "TE-(A)n", "TE-AluYk3",
        "TE-LTR7B", "TE-L1MCa", "TE-(GAGC)n", "TE-G-rich", "TE-Charlie1b")



for (x in te){
  p <- FeaturePlot(epit, features = x, reduction = "umap")
  file_name <- paste0("Epit_TEmarkers_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#DotPlot

features_to_plot <- c("TE-Harlequin-int", "TE-LTR48B", "TE-HERVH-int", "TE-(A)n", "TE-AluYk3",
                      "TE-LTR7B", "TE-L1MCa", "TE-(GAGC)n", "TE-G-rich", "TE-Charlie1b")

p <- DotPlot(object = epit, features = features_to_plot, cols = "RdBu", 
             group.by = "annotation", assay='TE') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("Epit_DotPlot_TEmarkers.pdf"), plot = combined_plot, width = 10, height = 8)

#FindMarkers (Basal vs Classical) -> no results

epit[["moffitt"]] <- NA
epit$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", epit$sample), "Basal",
                       ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", epit$sample), "Classical",
                              "Intermediate"))

epit.markers <- FindAllMarkers(epit, only.pos = TRUE, assay = "TE",
                               logfc.threshold = 1,min.pct=0.2, group.by = "moffitt")


top15_markers <- epit.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "Epit_TE_Top15_Markers_per_moffitt.csv", row.names = FALSE)




#-------------------------------     CAFs   ------------------------------------


Idents(seurat) <- "cell type"
DefaultAssay(seurat) <- "RNA"

#### Subset ####

caf <- subset(seurat, idents = "CAFs")

#...................................    RNA    .................................

#### Normalization and Scaling ####

#Normalize
caf <- NormalizeData(object = caf, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
caf <- FindVariableFeatures(object = caf, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
caf <- ScaleData(object = caf, features = VariableFeatures(caf))


#### Run SCT ####

caf <- SCTransform(caf, vars.to.regress = "percent.mt", verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

caf <- RunPCA(object = caf, features = VariableFeatures(object = caf))


#### Run Harmony and Determine Dimensions for 90% Variance ####

caf <- RunHarmony(caf, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- caf@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#### Find Neighbors e Clusters ####

#Find Neighbors 

caf <- FindNeighbors(caf, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

#Find Clusters

caf <- FindClusters(caf, verbose = TRUE, resolution = c(0, 0.05, 0.1, 0.2, 0.4, 0.6))

p <- clustree(caf)
ggsave("CAF_clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

caf <- RunUMAP(caf, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

saveRDS(caf,"Caf_UMAP.rds")----------------------------------------------------
  
Idents(caf) <- "SCT_snn_res.0.2"

#UMAP clusters

p <- DimPlot(caf, reduction = "umap")
ggsave("CAF_UMAP_clusters.pdf", plot = p , width = 10, height = 8)

#UMAP patients 

p <- DimPlot(caf, reduction = "umap", group.by = "sample")
ggsave("CAF_UMAP_samples.pdf", plot = p , width = 10, height = 8)

#UMAP Moffitt

caf[["moffitt"]] <- NA
caf$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", caf$sample), "Basal",
                       ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", caf$sample), "Classical",
                              "Intermediate"))


p <- DimPlot(caf, reduction = "umap", group.by = "SCT_snn_res.0.2", split.by = "moffitt")
ggsave("CAF_UMAP_moffitt.pdf", plot = p , width = 14, height = 8)


#### DotPlot ####

features_to_plot <- c("MMP11", "IGFBP3", "COL11A1", "COL10A1", "CTHRC1", "CFD", "PLA2G2A", "PTGDS", "C3", "C7",
                      "RGS5", "MYH11", "ADIRF", "MCAM", "MUSTN1", "CENPF", "H2AFZ", "STMN1", "TGFBI", "PTTG1",
                      "S100B", "CDH19", "GPM6B", "NRXN1", "CRYAB", "LYZ", "PGC", "TFF2", "LCN2", "KRT19",
                      "SPP1", "IBSP", "MMP13", "GZMA", "FABP5")

p <- DotPlot(object = caf, features = features_to_plot, cols = "RdBu", 
             group.by = "SCT_snn_res.0.2", assay='SCT') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("CAF_DotPlot.pdf"), plot = combined_plot, width = 10, height = 8)

#### FindMarkers ####

caf <- PrepSCTFindMarkers(caf)

caf.markers <- FindAllMarkers(caf, only.pos = TRUE, assay = "SCT",
                               logfc.threshold = 1,min.pct=0.2, group.by = "SCT_snn_res.0.2")

top15_markers <- caf.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "CAF_genes_Top15_Markers_per_Cluster.csv", row.names = FALSE)

#### COL3A1 VlnPlot ####

p <- VlnPlot(caf, features = "COL3A1")
ggsave("CAF_VlnPlot_COL3A1.pdf", plot = p , width = 10, height = 8)

#### Annotation Column ####

clusters <- FetchData(caf, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "myCAFs",
  "1" = "iCAFs",
  "2" = "Epithelial-like",
  "3" = "iCAFs",
  "4" = "Pericytes",
  "5" = "other CAFs",
  "6" = "CAFs (epigenetic state)",
  "7" = "Chondrocyte-like",
  "8" = "T/NK",
  "9" = "Erythrocytes",
  "10" = "CAFs (stressed)",
  "11" = "Peri-islet Schwann cells",
  "12" = "CAFs (high mitochondrial activity)",
  "13" = "Endothelial",
  "14" = "Epithelial"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

caf[["annotation"]] <- celltypes

saveRDS(caf,"CAF_celltypes.rds")-----------------------------------------------

  
#UMAP group.by cell type
  
Idents(caf) <- "SCT_snn_res.0.2"  
p <- DimPlot(caf, reduction = "umap", group.by = "annotation")
ggsave("CAF_UMAP_group.by_celltypes.pdf", plot = p , width = 10, height = 8)


#...................................    TEs    .................................


DefaultAssay(caf) <- "TE"


#### Normalization and Scaling ####

#Normalize
caf <- NormalizeData(object = caf, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
caf <- FindVariableFeatures(object = caf, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
caf <- ScaleData(object = caf, features = VariableFeatures(caf))

saveRDS(caf,"CAF_TE.rds")-----------------------------------------------------
  
  
#### TOP 15 TEs ####

te_matrix <- GetAssayData(caf, assay = "TE", layer = "data")

te_means <- rowMeans(as.matrix(te_matrix))

top15_tes <- sort(te_means, decreasing = TRUE)[1:15]

print(top15_tes)

#Feature Plot

top15 <- c("TE-AluSx", "TE-AluSz", "TE-AluSx1", "TE-AluY", "TE-AluJb", "TE-AluSp", 
           "TE-AluJr", "TE-L2a", "TE-AluSq2", "TE-AluSg", "TE-MIRb", "TE-AluSc", "TE-MIR",
           "TE-AluSx3", "TE-AluSz6")



for (x in top15){
  p <- FeaturePlot(caf, features = x, reduction = "umap")
  file_name <- paste0("CAF_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)

}

#### FindMarkers ####

#FindMarkers (annotation)

caf.markers <- FindAllMarkers(caf, only.pos = TRUE, assay = "TE",
                               logfc.threshold = 1,min.pct=0.2, group.by = "annotation")


top15_markers <- caf.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "CAF_TE_Top15_Markers_per_celltype.csv", row.names = FALSE)

# Feature Plot

te <- c("TE-L1PB1", "TE-(T)n", "TE-AluYm1", "TE-FRAM", "TE-AluSz6",
        "TE-LTR48B", "TE-(A)n", "TE-L1MA6", "TE-MER61-int", "TE-MLT1M")



for (x in te){
  p <- FeaturePlot(caf, features = x, reduction = "umap")
  file_name <- paste0("CAF_TEmarkers_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#DotPlot

features_to_plot <- c("TE-L1PB1", "TE-(T)n", "TE-AluYm1", "TE-FRAM", "TE-AluSz6",
                      "TE-LTR48B", "TE-(A)n", "TE-L1MA6", "TE-MER61-int", "TE-MLT1M")

p <- DotPlot(object = caf, features = features_to_plot, cols = "RdBu", 
             group.by = "annotation", assay='TE') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("CAF_DotPlot_TEmarkers.pdf"), plot = combined_plot, width = 10, height = 8)


#FindMarkers (Basal vs Classical) -> no results

caf[["moffitt"]] <- NA
caf$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", caf$sample), "Basal",
                       ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", caf$sample), "Classical",
                              "Intermediate"))

caf.markers <- FindAllMarkers(caf, only.pos = TRUE, assay = "TE",
                               logfc.threshold = 1,min.pct=0.2, group.by = "moffitt")


top15_markers <- caf.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "CAF_TE_Top15_Markers_per_moffitt.csv", row.names = FALSE)


#### CAFs wc (to understand if other cell types are skewing the results) ####

Idents(caf) <- "annotation"

caf_wc <- subset(caf, idents = c("myCAFs", "iCAFs", "Epithelial-like", "Pericytes", "other CAFs",
                                 "CAFs (epigenetic state)", "Chondrocyte-like", "CAFs (stressed)",
                                "Peri-islet Schwann cells", "CAFs (high mitochondrial activity)"))
#Normalization and scaling

#Normalize
caf_wc <- NormalizeData(object = caf_wc, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
caf_wc <- FindVariableFeatures(object = caf_wc, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
caf_wc <- ScaleData(object = caf_wc, features = VariableFeatures(caf_wc))

#DotPlot

features_to_plot <- c("TE-L1PB1", "TE-(T)n", "TE-AluYm1", "TE-FRAM", "TE-AluSz6",
                      "TE-LTR48B", "TE-(A)n", "TE-L1MA6", "TE-MER61-int", "TE-MLT1M")

p <- DotPlot(object = caf_wc, features = features_to_plot, cols = "RdBu", 
             group.by = "annotation", assay='TE') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("CAFWC_DotPlot_TEmarkers.pdf"), plot = combined_plot, width = 10, height = 8)


#---------------------------- CAFs wc ------------------------------------------

Idents(caf) <- "annotation"

caf_wc <- subset(caf, idents = c("myCAFs", "iCAFs", "Epithelial-like", "Pericytes", "other CAFs",
                                 "CAFs (epigenetic state)", "Chondrocyte-like", "CAFs (stressed)",
                                 "Peri-islet Schwann cells", "CAFs (high mitochondrial activity)"))

#...................................    RNA    .................................

#### Normalization and Scaling ####

#Normalize
caf_wc <- NormalizeData(object = caf_wc, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
caf_wc <- FindVariableFeatures(object = caf_wc, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
caf_wc <- ScaleData(object = caf_wc, features = VariableFeatures(caf_wc))


#### Run SCT ####

caf_wc <- SCTransform(caf_wc, vars.to.regress = "percent.mt", verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

caf_wc <- RunPCA(object = caf_wc, features = VariableFeatures(object = caf_wc))


#### Run Harmony and Determine Dimensions for 90% Variance ####

caf_wc <- RunHarmony(caf_wc, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- caf_wc@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#### Find Neighbors e Clusters ####

#Find Neighbors 

caf_wc <- FindNeighbors(caf_wc, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

#Find Clusters

caf_wc <- FindClusters(caf_wc, verbose = TRUE, resolution = c(0, 0.05, 0.1, 0.2, 0.4, 0.6))

p <- clustree(caf_wc)
ggsave("CAFWC_clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

caf_wc <- RunUMAP(caf_wc, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

saveRDS(caf_wc,"CafWC_UMAP.rds")----------------------------------------------------
  
Idents(caf_wc) <- "SCT_snn_res.0.2"

#UMAP clusters

p <- DimPlot(caf_wc, reduction = "umap")
ggsave("CAFWC_UMAP_clusters.pdf", plot = p , width = 10, height = 8)

#UMAP patients 

p <- DimPlot(caf_wc, reduction = "umap", group.by = "sample")
ggsave("CAFWC_UMAP_samples.pdf", plot = p , width = 10, height = 8)

#UMAP Moffitt

caf_wc[["moffitt"]] <- NA
caf_wc$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", caf$sample), "Basal",
                      ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", caf$sample), "Classical",
                             "Intermediate"))


p <- DimPlot(caf_wc, reduction = "umap", group.by = "SCT_snn_res.0.2", split.by = "moffitt")
ggsave("CAFWC_UMAP_moffitt.pdf", plot = p , width = 14, height = 8)


#### DotPlot ####

features_to_plot <- c("MMP11", "IGFBP3", "COL11A1", "COL10A1", "CTHRC1", "CFD", "PLA2G2A", "PTGDS", "C3", "C7",
                      "RGS5", "MYH11", "ADIRF", "MCAM", "MUSTN1", "CENPF", "H2AFZ", "STMN1", "TGFBI", "PTTG1",
                      "S100B", "CDH19", "GPM6B", "NRXN1", "CRYAB", "LYZ", "PGC", "TFF2", "LCN2", "KRT19",
                      "SPP1", "IBSP", "MMP13", "GZMA", "FABP5")

p <- DotPlot(object = caf_wc, features = features_to_plot, cols = "RdBu", 
             group.by = "SCT_snn_res.0.2", assay='SCT') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("CAFWC_DotPlot.pdf"), plot = combined_plot, width = 10, height = 8)

#### FindMarkers ####

caf_wc <- PrepSCTFindMarkers(caf_wc)

caf_wc.markers <- FindAllMarkers(caf_wc, only.pos = TRUE, assay = "SCT",
                              logfc.threshold = 1,min.pct=0.2, group.by = "SCT_snn_res.0.2")

top15_markers <- caf_wc.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "CAFWC_genes_Top15_Markers_per_Cluster.csv", row.names = FALSE)

#### COL3A1 VlnPlot ####

p <- VlnPlot(caf_wc, features = "COL3A1")
ggsave("CAFWC_VlnPlot_COL3A1.pdf", plot = p , width = 10, height = 8)

#### Annotation Column ####

clusters <- FetchData(caf_wc, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "myCAFs",
  "1" = "Other CAFs",
  "2" = "Epithelial-like",
  "3" = "iCAFs",
  "4" = "Pericytes",
  "5" = "iCAFs",
  "6" = "CAFs (epigenetic state)",
  "7" = "Chondrocyte-like",
  "8" = "iCAFs",
  "9" = "Peri-islet Schwann cells",
  "10" = "iCAFs",
  "11" = "Endothelial"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

caf_wc[["annotation"]] <- celltypes

saveRDS(caf_wc,"CAF_celltypes.rds")-----------------------------------------------
  
#### BarPlot ####

#subset and column

caf_type <- subset(caf_wc, subset = SCT_snn_res.0.2 %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

clusters <- FetchData(caf_type, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "myCAFs",
  "1" = "Other CAFs",
  "2" = "Epithelial-like",
  "3" = "iCAFs",
  "4" = "Pericytes",
  "5" = "iCAFs",
  "6" = "CAFs (epigenetic state)",
  "7" = "Chondrocyte-like",
  "8" = "iCAFs",
  "9" = "Peri-islet Schwann cells",
  "10" = "iCAFs"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

caf_type[["CAF types"]] <- celltypes

#barplot

plot_data <- caf_type@meta.data %>%
  group_by(orig.ident, `CAF types`) %>% 
  summarise(cnt = n()) %>%
  group_by(orig.ident) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  ungroup()

p <- ggplot(plot_data, aes(x = orig.ident, y = freq, fill = `CAF types`)) +
  geom_bar(stat = "identity", width = 0.8) + 
  theme_bw() + 
  labs(
    y = "Proportion of cells in cluster", 
    x = NULL, 
    fill = NULL
  ) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "bottom", 
    panel.grid = element_blank(), 
    axis.line = element_line(colour = "black")
  )

ggsave("CAF_barplot.pdf", plot = p , width = 10, height = 8)#barplot  
#UMAP group.by cell type
  
Idents(caf_wc) <- "SCT_snn_res.0.2"  
p <- DimPlot(caf_wc, reduction = "umap", group.by = "annotation")
ggsave("CAFWC_UMAP_group.by_celltypes.pdf", plot = p , width = 10, height = 8)


#...................................    TEs    .................................


DefaultAssay(caf_wc) <- "TE"


#### Normalization and Scaling ####

#Normalize
caf_wc <- NormalizeData(object = caf_wc, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
caf_wc <- FindVariableFeatures(object = caf_wc, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
caf_wc <- ScaleData(object = caf_wc, features = VariableFeatures(caf_wc))

saveRDS(caf_wc,"CAF_TE.rds")-----------------------------------------------------
  
  
#### TOP 15 TEs ####

te_matrix <- GetAssayData(caf_wc, assay = "TE", layer = "data")

te_means <- rowMeans(as.matrix(te_matrix))

top15_tes <- sort(te_means, decreasing = TRUE)[1:15]

print(top15_tes)

#Feature Plot

top15 <- c("TE-AluSx", "TE-AluSz", "TE-AluSx1", "TE-AluY", "TE-AluJb", "TE-AluSp", 
           "TE-AluJr", "TE-L2a", "TE-AluSq2", "TE-AluSg", "TE-MIRb", "TE-AluSc", "TE-MIR",
           "TE-AluSx3", "TE-AluSz6")



for (x in top15){
  p <- FeaturePlot(caf_wc, features = x, reduction = "umap")
  file_name <- paste0("CAFWC_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#### FindMarkers ####

#FindMarkers (annotation)

caf_wc.markers <- FindAllMarkers(caf_wc, only.pos = TRUE, assay = "TE",
                              logfc.threshold = 1,min.pct=0.2, group.by = "annotation")


top15_markers <- caf_wc.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "CAFWC_TE_Top15_Markers_per_celltype.csv", row.names = FALSE)

# Feature Plot

te <- c("TE-L1PB1")



for (x in te){
  p <- FeaturePlot(caf_wc, features = x, reduction = "umap", split.by = "moffitt")
  file_name <- paste0("CAFWC_TEmarkers_FeaturePlot_moffitt_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#DotPlot

features_to_plot <- c("TE-L1PB1")

p <- DotPlot(object = caf_wc, features = features_to_plot, cols = "RdBu", 
             group.by = "annotation", assay='TE') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("CAFWC_DotPlot_TEmarkers.pdf"), plot = combined_plot, width = 10, height = 8)


#FindMarkers (Basal vs Classical) -> no results

caf_wc[["moffitt"]] <- NA
caf_wc$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", caf_wc$sample), "Basal",
                      ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", caf_wc$sample), "Classical",
                             "Intermediate"))

caf_wc.markers <- FindAllMarkers(caf_wc, only.pos = TRUE, assay = "TE",
                              logfc.threshold = 1,min.pct=0.2, group.by = "moffitt")


top15_markers <- caf_wc.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "CAFWC_TE_Top15_Markers_per_moffitt.csv", row.names = FALSE)


#### COL3A1 VlnPlot ####

p <- VlnPlot(caf_wc, features = "TE-L1PB1", split.by = "moffitt")
ggsave("CAFWC_VlnPlot_TE-L1PB1_moffitt.pdf", plot = p , width = 10, height = 8)




#---------------------------- Myeloid ------------------------------------------

Idents(seurat) <- "cell type"
DefaultAssay(seurat) <- "RNA"

#### Subset ####
mye <- subset(seurat, idents = "Myeloid")

#...................................    RNA    .................................

#### Normalization and Scaling ####

#Normalize
mye <- NormalizeData(object = mye, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
mye <- FindVariableFeatures(object = mye, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
mye <- ScaleData(object = mye, features = VariableFeatures(mye))


#### Run SCT ####

mye <- SCTransform(mye, vars.to.regress = "percent.mt", verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

mye <- RunPCA(object = mye, features = VariableFeatures(object = mye))


#### Run Harmony and Determine Dimensions for 90% Variance ####

mye <- RunHarmony(mye, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- mye@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#### Find Neighbors e Clusters ####

#Find Neighbors 

mye <- FindNeighbors(mye, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

#Find Clusters

mye <- FindClusters(mye, verbose = TRUE, resolution = c(0, 0.05, 0.1, 0.2, 0.4, 0.6))

p <- clustree(mye)
ggsave("MYE_clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

mye <- RunUMAP(mye, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

saveRDS(mye,"MYE_UMAP.rds")----------------------------------------------------
  
Idents(mye) <- "SCT_snn_res.0.2"

#UMAP clusters

p <- DimPlot(mye, reduction = "umap")
ggsave("MYE_UMAP_clusters.pdf", plot = p , width = 10, height = 8)

#UMAP patients 

p <- DimPlot(mye, reduction = "umap", group.by = "sample")
ggsave("MYE_UMAP_samples.pdf", plot = p , width = 10, height = 8)

#UMAP Moffitt

mye[["moffitt"]] <- NA
mye$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", mye$sample), "Basal",
                      ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", mye$sample), "Classical",
                             "Intermediate"))


p <- DimPlot(mye, reduction = "umap", group.by = "SCT_snn_res.0.2", split.by = "moffitt")
ggsave("MYE_UMAP_moffitt.pdf", plot = p , width = 14, height = 8)


#### DotPlot ####

features_to_plot <- c("SPP1", "MARCO", "C1QA", "C1QB", "C1QC", "S100A8", "S100A9", "S100A12",
                      "FCGR3A", "CDKN1C", "CLEC9A", "XCR1", "CD1C", "FCER1A", "CCL19", "CCL22",
                      "CCR7", "LILRA4", "PLD4", "MKI67", "KIT", "CD3E", "KRT18", "KRT19")

p <- DotPlot(object = mye, features = features_to_plot, cols = "RdBu", 
             group.by = "SCT_snn_res.0.2", assay='SCT') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("MYE_DotPlot.pdf"), plot = combined_plot, width = 10, height = 8)

#### FindMarkers ####

mye <- PrepSCTFindMarkers(mye)

mye.markers <- FindAllMarkers(mye, only.pos = TRUE, assay = "SCT",
                              logfc.threshold = 1,min.pct=0.2, group.by = "SCT_snn_res.0.2")

top15_markers <- mye.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "MYE_genes_Top15_Markers_per_Cluster.csv", row.names = FALSE)

#### CD68 VlnPlot ####

p <- VlnPlot(mye, features = "CD68")
ggsave("MYE_VlnPlot_CD68.pdf", plot = p , width = 10, height = 8)

#### Annotation Column ####

clusters <- FetchData(mye, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "C1QC+ macrophage",
  "1" = "SPP1+ macrophage",
  "2" = "MDSC",
  "3" = "MDSC",
  "4" = "Epithelial-like",
  "5" = "cDC2",
  "6" = "Myeloid (inflamatory)",
  "7" = "Myeloid (stressed)",
  "8" = "T",
  "9" = "Erythrocytes",
  "10" = "Monocyte",
  "11" = "CAFS",
  "12" = "Proliferating",
  "13" = "Myeloid (stressed)",
  "14" = "Epithelial"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

mye[["annotation"]] <- celltypes

saveRDS(mye,"MYE_celltypes.rds")-----------------------------------------------
  
  
#UMAP group.by cell type
  
Idents(mye) <- "SCT_snn_res.0.2"  
p <- DimPlot(mye, reduction = "umap", group.by = "annotation")
ggsave("MYE_UMAP_group.by_celltypes.pdf", plot = p , width = 10, height = 8)


#...................................    TEs    .................................


DefaultAssay(mye) <- "TE"


#### Normalization and Scaling ####

#Normalize
mye <- NormalizeData(object = mye, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
mye <- FindVariableFeatures(object = mye, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
mye <- ScaleData(object = mye, features = VariableFeatures(mye))

saveRDS(mye,"MYE_TE.rds")-----------------------------------------------------
  
  
#### TOP 15 TEs ####

te_matrix <- GetAssayData(mye, assay = "TE", layer = "data")

te_means <- rowMeans(as.matrix(te_matrix))

top15_tes <- sort(te_means, decreasing = TRUE)[1:15]

print(top15_tes)

#Feature Plot

top15 <- c("TE-AluSx", "TE-AluSz", "TE-AluSx1", "TE-AluY", "TE-AluJb", "TE-AluSp", 
           "TE-AluJr", "TE-L2a", "TE-AluSq2", "TE-AluSg", "TE-MIRb", "TE-AluSc", "TE-MIR",
           "TE-AluSx3", "TE-AluSz6")



for (x in top15){
  p <- FeaturePlot(mye, features = x, reduction = "umap")
  file_name <- paste0("MYE_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#### FindMarkers ####

#FindMarkers (annotation)

mye.markers <- FindAllMarkers(mye, only.pos = TRUE, assay = "TE",
                              logfc.threshold = 1,min.pct=0.2, group.by = "annotation")


top15_markers <- mye.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "MYE_TE_Top15_Markers_per_celltype.csv", row.names = FALSE)

# Feature Plot

te <- c("TE-(T)n", "TE-MLT1D", "TE-L1MB7", "TE-AluSg4", "TE-L1MEd",
        "TE-G-rich", "TE-MER21C", "TE-L1PA8", "TE-L1MCa", "TE-L1M2")



for (x in te){
  p <- FeaturePlot(mye, features = x, reduction = "umap")
  file_name <- paste0("MYE_TEmarkers_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#DotPlot

features_to_plot <- c("TE-(T)n", "TE-MLT1D", "TE-L1MB7", "TE-AluSg4", "TE-L1MEd",
                      "TE-G-rich", "TE-MER21C", "TE-L1PA8", "TE-L1MCa", "TE-L1M2")

p <- DotPlot(object = mye, features = features_to_plot, cols = "RdBu", 
             group.by = "annotation", assay='TE') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("MYE_DotPlot_TEmarkers.pdf"), plot = combined_plot, width = 10, height = 8)


#-------------------------- Myeloid wc -----------------------------------------

Idents(mye) <- "annotation"

mye <- subset(mye, idents = c("C1QC+ macrophage", "SPP1+ macrophage", "MDSC", "Epithelial-like", "cDC2",
                                 "Myeloid (inflamatory)", "Myeloid (stressed)", "Monocyte",
                                 "Proliferating"))

#...................................    RNA    .................................

#### Normalization and Scaling ####

#Normalize
mye <- NormalizeData(object = mye, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
mye <- FindVariableFeatures(object = mye, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
mye <- ScaleData(object = mye, features = VariableFeatures(mye))


#### Run SCT ####

mye <- SCTransform(mye, vars.to.regress = "percent.mt", verbose = TRUE, conserve.memory = TRUE)

#Run PCA 

mye <- RunPCA(object = mye, features = VariableFeatures(object = mye))


#### Run Harmony and Determine Dimensions for 90% Variance ####

mye <- RunHarmony(mye, group.by.vars = c("sample"), plot_convergence = F) 

stdev <- mye@reductions$harmony@stdev 
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}

#Confirm PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#### Find Neighbors e Clusters ####

#Find Neighbors 

mye <- FindNeighbors(mye, dims = 1:PCNum, verbose = TRUE, reduction = "harmony")

#Find Clusters

mye <- FindClusters(mye, verbose = TRUE, resolution = c(0, 0.05, 0.1, 0.2, 0.4, 0.6))

p <- clustree(mye)
ggsave("MYEWC_clustree.pdf", plot = p , width = 10, height = 8)


#### UMAPs ####

mye <- RunUMAP(mye, dims = 1:PCNum, verbose = TRUE, reduction = "harmony", n.components = 2L)

saveRDS(mye,"MYEWC_UMAP.rds")----------------------------------------------------
  
Idents(mye) <- "SCT_snn_res.0.2"

#UMAP clusters

p <- DimPlot(mye, reduction = "umap")
ggsave("MYEWC_UMAP_clusters.pdf", plot = p , width = 10, height = 8)

#UMAP patients 

p <- DimPlot(mye, reduction = "umap", group.by = "sample")
ggsave("MYEWC_UMAP_samples.pdf", plot = p , width = 10, height = 8)

#UMAP Moffitt

mye[["moffitt"]] <- NA
mye$moffitt <- ifelse(grepl("^(P11|P14|P17|P18)", mye$sample), "Basal",
                      ifelse(grepl("^(P01|P04|P05|P07|P08|P09|P12|P16|P20|P23)", mye$sample), "Classical",
                             "Intermediate"))


p <- DimPlot(mye, reduction = "umap", group.by = "SCT_snn_res.0.2", split.by = "moffitt")
ggsave("MYEWC_UMAP_moffitt.pdf", plot = p , width = 14, height = 8)


#### DotPlot ####

features_to_plot <- c("SPP1", "MARCO", "C1QA", "C1QB", "C1QC", "S100A8", "S100A9", "S100A12",
                      "FCGR3A", "CDKN1C", "CLEC9A", "XCR1", "CD1C", "FCER1A", "CCL19", "CCL22",
                      "CCR7", "LILRA4", "PLD4", "MKI67", "KIT", "CD3E", "KRT18", "KRT19")

p <- DotPlot(object = mye, features = features_to_plot, cols = "RdBu", 
             group.by = "SCT_snn_res.0.2", assay='SCT') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("MYEWC_DotPlot.pdf"), plot = combined_plot, width = 10, height = 8)

#### FindMarkers ####

mye <- PrepSCTFindMarkers(mye)

mye.markers <- FindAllMarkers(mye, only.pos = TRUE, assay = "SCT",
                              logfc.threshold = 1,min.pct=0.2, group.by = "SCT_snn_res.0.2")

top15_markers <- mye.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "MYEWC_genes_Top15_Markers_per_Cluster.csv", row.names = FALSE)

#### CD68 VlnPlot ####

p <- VlnPlot(mye, features = "CD68")
ggsave("MYEWC_VlnPlot_CD68.pdf", plot = p , width = 10, height = 8)

#### Annotation Column ####

clusters <- FetchData(mye, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "cDC2",
  "1" = "C1QC+ macrophage",
  "2" = "SPP1+ macrophage",
  "3" = "MDSC",
  "4" = "MDSC",
  "5" = "Epithelial-like",
  "6" = "Myeloid (inflamatory)",
  "7" = "Myeloid (stressed)",
  "8" = "Monocyte",
  "9" = "C1QC+ macrophage",
  "10" = "Proliferating",
  "11" = "Myeloid (stressed)",
  "12" = "cDC3"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

mye[["annotation"]] <- celltypes

saveRDS(mye,"MYEWC_celltypes.rds")-----------------------------------------------
  
  
#UMAP group.by cell type
  
Idents(mye) <- "SCT_snn_res.0.2"  
p <- DimPlot(mye, reduction = "umap", group.by = "annotation")
ggsave("MYEWC_UMAP_group.by_celltypes.pdf", plot = p , width = 10, height = 8)

#### BarPlot ####

#subset and column

mye_type <- subset(mye, subset = SCT_snn_res.0.2 %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

clusters <- FetchData(mye_type, vars = "SCT_snn_res.0.2")

cell_map <- c(
  "0" = "cDC2",
  "1" = "C1QC+ macrophage",
  "2" = "SPP1+ macrophage",
  "3" = "MDSC",
  "4" = "MDSC",
  "5" = "Epithelial-like",
  "6" = "Myeloid (inflamatory)",
  "7" = "Myeloid (stressed)",
  "8" = "Monocyte",
  "9" = "C1QC+ macrophage",
  "10" = "Proliferating",
  "11" = "Myeloid (stressed)",
  "12" = "cDC3"
)

celltypes <- cell_map[as.character(clusters$SCT_snn_res.0.2)] 

names(celltypes) <- rownames(clusters)

mye_type[["MYE types"]] <- celltypes

#barplot

plot_data <- mye_type@meta.data %>%
  group_by(orig.ident, `MYE types`) %>% 
  summarise(cnt = n()) %>%
  group_by(orig.ident) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  ungroup()

p <- ggplot(plot_data, aes(x = orig.ident, y = freq, fill = `MYE types`)) +
  geom_bar(stat = "identity", width = 0.8) + 
  theme_bw() + 
  labs(
    y = "Proportion of cells in cluster", 
    x = NULL, 
    fill = NULL
  ) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
    legend.position = "bottom", 
    panel.grid = element_blank(), 
    axis.line = element_line(colour = "black")
  )

ggsave("MYEWC_barplot.pdf", plot = p , width = 10, height = 8)#barplot  
#...................................    TEs    .................................


DefaultAssay(mye) <- "TE"


#### Normalization and Scaling ####

#Normalize
mye <- NormalizeData(object = mye, normalization.method = "LogNormalize", scale.factor = 10000)

#Variable Features
mye <- FindVariableFeatures(object = mye, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures = 5000)

#Scale Data
mye <- ScaleData(object = mye, features = VariableFeatures(mye))

saveRDS(mye,"MYEWC_TE.rds")-----------------------------------------------------
  
  
#### TOP 15 TEs ####

te_matrix <- GetAssayData(mye, assay = "TE", layer = "data")

te_means <- rowMeans(as.matrix(te_matrix))

top15_tes <- sort(te_means, decreasing = TRUE)[1:15]

print(top15_tes)

#Feature Plot

top15 <- c("TE-AluSx", "TE-AluSz", "TE-AluSx1", "TE-AluY", "TE-AluJb", "TE-AluSp", 
           "TE-AluJr", "TE-L2a", "TE-AluSq2", "TE-AluSg", "TE-MIRb", "TE-AluSc", "TE-MIR",
           "TE-AluSx3", "TE-AluSz6")



for (x in top15){
  p <- FeaturePlot(mye, features = x, reduction = "umap")
  file_name <- paste0("MYEWC_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#### FindMarkers ####

#FindMarkers (annotation)

mye.markers <- FindAllMarkers(mye, only.pos = TRUE, assay = "TE",
                              logfc.threshold = 1,min.pct=0.2, group.by = "annotation")


top15_markers <- mye.markers %>% group_by(cluster) %>% slice_max(n = 15, order_by = avg_log2FC) 
write.csv(top15_markers, "MYEWC_TE_Top15_Markers_per_celltype.csv", row.names = FALSE)

# Feature Plot

te <- c("TE-(T)n", "TE-MLT1D", "TE-L1MB7", "TE-L1MEd", "TE-MLT1D",
        "TE-HERVH-int", "TE-G-rich", "TE-L1MA3", "TE-L1PA8", "TE-L1M2", "TE-L1MCa",
        "TE-MER34B-int", "TE-L1M2", "TE-L1ME3D", "TE-L1MA4", "TE-MLT1H")



for (x in te){
  p <- FeaturePlot(mye, features = x, reduction = "umap")
  file_name <- paste0("MYEWC_TEmarkers_FeaturePlot_", x, ".pdf")
  ggsave(file_name, plot = p , width = 10, height = 8)
  
}

#DotPlot

features_to_plot <- c("TE-(T)n", "TE-MLT1D", "TE-L1MB7", "TE-L1MEd",
                      "TE-HERVH-int", "TE-G-rich", "TE-L1MA3", "TE-L1PA8", "TE-L1M2", "TE-L1MCa",
                      "TE-MER34B-int", "TE-L1ME3D", "TE-L1MA4", "TE-MLT1H")

p <- DotPlot(object = mye, features = features_to_plot, cols = "RdBu", 
             group.by = "annotation", assay='TE') + RotatedAxis()

combined_plot <- wrap_plots(p)

ggsave(paste0("MYEWC_DotPlot_TEmarkers.pdf"), plot = combined_plot, width = 10, height = 8)



# scRNA-seq lab
# BCB 5250
# Megan Johnson
# Clustering

# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
seurat_integrated <- readRDS("results/integrated_seurat.rds")

## Identifying significant PCs:

# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)


## Clustering the cells:

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


## Visualizing clusters of cells

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

# For resolution of 0.8:
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# For resolution of 0.4:
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# My clusters look different from those in the tutorial, so I am loading their result
# before proceeding to the next steps:
load(bzfile("data/additional_data/seurat_integrated.RData.bz2"))

## Checking the object clusters with different resolutions:

# For resolution of 0.4:
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# There are 13 clusters at resolution 0.4. This may not be enough - it looks like some of the
# clusters could be subdivided.

# For resolution of 0.6:
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# There are 15 clusters at resolution 0.6. This looks pretty good.

# For resolution of 0.8:
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# There are 17 clusters at resolution 0.8. This also looks pretty good. The clusters
# make sense.

# For resolution of 1.0:
Idents(object = seurat_integrated) <- "integrated_snn_res.1"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# There are 22 clusters at resolution 1.0. This may be too many clusters. We can see that some of
# them only represent a very small number of points.

# For resolution of 1.4:
Idents(object = seurat_integrated) <- "integrated_snn_res.1.4"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# There are 27 clusters at resolution 1.4. This looks like too many clusters. 


## Proceeding with resolution 0.8.
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


## Clustering Quality Control

# QC Metrics:

# Segregation of clusters by sample
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)
# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), position=position_fill())
# These visualizations show that the clusters are similar between the ctrl and stim conditions,
# which is what we expect.

# Segregation of clusters by cell cycle phase
# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
# The clustering looks very similar for each cell cycle phase.

# Segregation of clusters by other sources of uninteresting variation
# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, fill=integrated_snn_res.0.8)) +
  NoLegend()

# Exploration of the PCs driving the clusters
# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")
# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

# Exploration of known cell type markers
# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"
# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

# CD14+ monocyte markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# CD14+ monocyte markers are highly expressed in clusters 1 and 3.

# FCGR3A+ monocyte markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# FCGR3A+ monocyte markers are highly expressed in cluster 10.

# Macrophages
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Macrophage markers are not highly expressed in any specific cluster.

# Conventional dendritic cell markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Conventional dendritic cell markers are highly expressed in cluster 14.

# Plasmacytoid dendritic cell markers
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Plasmacytoid dendritic cell markers are highly expressed in cluster 16.

# Can also express expression across clusters using Seurat DotPlot visualization tool
# List of known celltype markers
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")
# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")

## Exercise: Hypothesizing the clusters corresponding to the remaining cell types.
markers[["B cells"]] <- c("CD79A", "MS4A1")
markers[["T cells"]] <- c("CD3D")
markers[["CD4+ T cells"]] <- c("CD3D", "IL7R", "CCR7")
markers[["CD8+ T cells"]] <- c("CD3D", "ICD8A")
markers[["NK cells"]] <- c("GNLY", "NKG7")
markers[["Megakaryocytes"]] <- c("PPBP")
markers[["Erythrocytes"]] <- c("HBB", "HBA2")

# B cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["B cells"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Clusters 7 and 11 correspond to B cells.

# T cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["T cells"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# This one is hard to say because there are multiple clusters where the T cell marker is highly
# expressed (clusters 5, 9, 2, 0, and 6). Having more than one marker gene would be helpful for
# resolving this.

# CD4+ T cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["CD4+ T cells"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Clusters 0, 2, 4, and 6 correspond to CD4+ T cells.

# CD8+ T cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["CD8+ T cells"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Cluster 5 appears to correspond with CD8+ T cells.

# NK cells
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["NK cells"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Clusters 8 and 12 correspond to NK cells.

# Megakaryocytes
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["Megakaryocytes"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# Cluster 15 corresponds to megakaryocytes.

# Erythrocytes
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = markers[["Erythrocytes"]], 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
# None of the clusters correspond to erythrocytes.

## Questions:
# 1. T cell markers appear to be highly expressed in many clusters. How can we differentiate and
# subset the larger group into smaller subset of cells?
# Ans: T cell markers are associated with clusters 5, 9, 2, 0, and 6. We may have to use more PCs
# to be able to differentiate between these.
# 2. Do the clusters corresponding to the same cell types have biologically meaningful differences?
# Are there subpopulations of these cell types?
# Ans: For some clusters corresponding to the same cell type, it seems likely that there are
# biologically meaningful differences. For example, clusters 7 and 11 both correspond to B cells
# but are clearly separate clusters on the UMAP. Clusters, 8 and 12, on the other hand, both
# correspond to NK cells but 8 and 12 could just be one big cluster on the UMAP.
# 3. Can we acquire higher confidence in these cell type identities by identifying other marker
# genes for these clusters?
# Ans: Yes, identifying other marker genes for the clusters would definitely help increase our
# confidence in the cell type identification.


## Marker identification

# Identification of conserved markers
DefaultAssay(seurat_integrated) <- "RNA"
# Testing on one cluster:
cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
# Adding gene annotation
annotations <- read.csv("data/annotation.csv")
# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
View(cluster0_ann_markers)

# Exercise: In the previous lesson, we identified cluster 10 as FCGR3A+ monocytes by inspecting
# the expression of known cell markers FCGR3A and MS4A7. Use FindConservedMarkers() function to
# find conserved markers for cluster 10. What do you observe? Do you see FCGR3A and MS4A7 as highly
# expressed genes in cluster 10?
cluster10_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 10,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
# FCGR3A and MS4A7 are both found to be highly expressed in cluster 10. MS4A4A, CXCL16, and VMO1
# were also found to be highly expressed in this cluster.

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
# Iterate function across CD4+ T cell clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_log2FC + stim_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)
# Visualize top 10 markers per cluster
View(top10)
# Plot interesting marker gene expression for cluster 4
FeaturePlot(object = seurat_integrated, 
            features = c("HSPH1", "HSPE1", "DNAJB1"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
# Vln plot - cluster 4
VlnPlot(object = seurat_integrated, 
        features = c("HSPH1", "HSPE1", "DNAJB1"))

# Identify genes that are differentially expressed between the markers
# Determine differentiating markers for CD4+ T cell
cd4_tcells <- FindMarkers(seurat_integrated,
                          ident.1 = 2,
                          ident.2 = c(0,4,6))                  
# Add gene symbols to the DE table
cd4_tcells <- cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))
# Reorder columns and sort by padj      
cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]
cd4_tcells <- cd4_tcells %>%
  dplyr::arrange(p_val_adj) 
# View data
View(cd4_tcells)


# Having identified cell types, we now can create a labeled UMAP
# Rename all identities
seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "Naive or memory CD4+ T cells",
                                  "1" = "CD14+ monocytes",
                                  "2" = "Activated T cells",
                                  "3" = "CD14+ monocytes",
                                  "4" = "Stressed cells / Unknown",
                                  "5" = "CD8+ T cells",
                                  "6" = "Naive or memory CD4+ T cells",
                                  "7" = "B cells",
                                  "8" = "NK cells",
                                  "9" = "CD8+ T cells",
                                  "10" = "FCGR3A+ monocytes",
                                  "11" = "B cells",
                                  "12" = "NK cells",
                                  "13" = "B cells",
                                  "14" = "Conventional dendritic cells",
                                  "15" = "Megakaryocytes",
                                  "16" = "Plasmacytoid dendritic cells")
# Plot the UMAP
DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
# Remove the stressed or dying cells
seurat_subset_labeled <- subset(seurat_integrated,
                                idents = "Stressed cells / Unknown", invert = TRUE)
# Re-visualize the clusters
DimPlot(object = seurat_subset_labeled, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)
# Save final labeled R object
write_rds(seurat_integrated,
          file = "results/seurat_labelled.rds")
# Create and save a text file with sessionInfo
sink("sessionInfo_scrnaseq_Feb2023.txt")
sessionInfo()
sink()

# Determine if there is a shift in cell populations between ctrl and stim.
# Add celltype annotation as a column in meta.data 
seurat_subset_labeled$celltype <- Idents(seurat_subset_labeled)
# Compute number of cells per celltype
n_cells <- FetchData(seurat_subset_labeled, 
                     vars = c("celltype", "sample")) %>%
  dplyr::count(celltype, sample)
# Barplot of number of cells per celltype by sample
ggplot(n_cells, aes(x=celltype, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))

# Perform differential expression analysis between conditions ctrl and stim
# Subset seurat object to just B cells
seurat_b_cells <- subset(seurat_subset_labeled, subset = (celltype == "B cells"))

# Run a wilcox test to compare ctrl vs stim
b_markers <- FindMarkers(seurat_b_cells@meta.data$sample,
                         ident.1 = "ctrl",
                         ident.2 = "stim",
                         only.pos = FALSE,
                         logfc.threshold = 0.25)

# Visualize using an enhanced volcano plot
library(EnhancedVolcano)
EnhancedVolcano(b_markers,
                row.names(b_markers),
                x="avg_log2FC",
                y="p_val_adj"
)

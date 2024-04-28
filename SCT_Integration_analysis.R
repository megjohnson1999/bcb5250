# scRNA-seq lab
# BCB 5250
# Megan Johnson
# Normalization and Integration

# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

# Evaluation for cell cycle phase:
# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)  

# Use PCA to evaluate the effects of cell cycle
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
# The PCA plots for each cell cycle phase look similar, so for our data we would not regress
# out the variation from cell cycle.


## Evaluation for mitochondrial expression
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))

# Plot the PCA colored by mitoRatio
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")
# The PCA plots for each mitoRatio category look similar, so I think we would not regress out the
# variation from mitoRatio.

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# Make sure R can handle large object size
options(future.globals.maxSize = 4000 * 1024^2)

# Loop to perform the sctransform on all samples
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
}

# Check which assays are stored in objects
split_seurat$ctrl@assays
split_seurat$stim@assays
# Question 1: Are the same assays available for the “stim” samples within the split_seurat object?
# What is the code you used to check that?
# Answer: See code above. The same assays (RNA and SCT) are available for both the ctrl and the
# stim samples.
# Question 2: Any observations for the genes or features listed under “First 10 features:” and the
# “Top 10 variable features:” for “ctrl” versus “stim”?
# Answer: Many genes appear in the top 10 variable features for both ctrl and stim.

# Save the split seurat object for future use
saveRDS(split_seurat, "data/split_seurat.rds")


## Integration

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Identify anchors
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Visualizing the integrated data using PCA
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)
# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample")  

# Visualizing the integrated data using UMAP
# Set seed
set.seed(123456)
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated)   
# Plot UMAP split by sample
DimPlot(seurat_integrated,
        split.by = "sample")  

# Save integrated seurat object
saveRDS(seurat_integrated, "results/integrated_seurat.rds")
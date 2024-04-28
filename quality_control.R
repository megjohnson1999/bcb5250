# scRNA-seq lab
# BCB 5250
# Megan Johnson
# Quality Control

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(SeuratWrappers)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)


# Single-cell RNA-seq analysis - QC

# Read in the data
# Create a Seurat object for each sample 
for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = file)
  assign(file, seurat_obj) }
# Question: What is min.features? Describe it and min.features = 100. 
# Answer: min.features is the minimum number of genes that need to be detected for a cell to be
# included. Here we use min.features = 100, so any cell with fewer than 100 genes detected is
# filtered out of the data.

# Check the metadata in the new Seurat objects 
head(ctrl_raw_feature_bc_matrix@meta.data) 
head(stim_raw_feature_bc_matrix@meta.data)
# Question: Report nCount_RNA and nFeature_RNA for AAACATACATTTCC-1 (2nd item) from the
# ctrl_raw_feature_bc_matrix.
# Answer: For AAACATACATTTCC-1, nCount is 3125 and nFeature is 896.

# Merge into a single Seurat object
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, 
                       y = stim_raw_feature_bc_matrix, 
                       add.cell.id = c("ctrl", "stim"))
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# Add number of genes per UMI for each cell to metadata (novelty score)
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Calculate mitochondrial ratio (proportion of reads originating from mitochondrial genes)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Extract the metadata into a new dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
# Rename columns to be more intuitive
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")


# Assessing the quality metrics

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# This bar graph shows over 15,000 cells per sample, which is more than the 12,000-13,000 that
# were expected for this experiment.

# Visualize the number of UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
# The vertical line shows nUMI=500, which is the lower cutoff for cells to be included.
# Cells with fewer than 500 transcripts should not be included in the analysis.

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# The vertical line shows nGene=300, which is the low end of what we are looking for for genes
# detected per cell.

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
# The vertical line shows novelty score=0.8. Good quality cells are expected to have novelty score
# of at least 0.8, and we see that this data looks very good for this metric.

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# The vertical line shows mitochondrial ratio=0.2. Mitochondrial ratio > 0.2 suggests contamination
# from dead or dying cells. We see that this data looks very good for this metric.

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)


## Filtering

# Filter out low quality cells using selected thresholds
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

# Identifying genes with zero counts
filtered_seurat <- JoinLayers(filtered_seurat)
counts <- LayerData(object=filtered_seurat, layer="counts")
# I used https://github.com/satijalab/seurat/issues/7905 to figure out what to use instead of
# GetAssayData(), which was not working because of issues with the layers.
nonzero <- counts > 0
# Filter out genes expressed in 10 or fewer cells
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

## QC metrics again with filtered data
# Visualize the number of cell counts per sample
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# After filtering we now have slightly less than 15000 cells per sample.

# Visualize the number UMIs/transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500, lty=2)
# Cells with nUMI < 500 were filtered out.

# Visualize the distribution of genes detected per cell via histogram
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 250, lty=2)
# Cells with nGene < 250 were filtered out.

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8, lty=2)
# Cells with novelty score <= 0.8 were filtered out.

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# Cells with mitoRatio >= 0.2 were filtered out.

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

## Questions
# Report the number of cells left for each sample, and comment on whether the number of cells
# removed is high or low. Can you give reasons why this number is still not ~12K (which is how
# many cells were loaded for the experiment)?
metadata_clean %>%
  count(sample)
# Ans: There are 14847 ctrl cells and 14782 stim cells remaining after filtering. The number of
# cells removed is low (more cells should have been removed to get the counts down to the
# expected 12000 per sample). This may mean that our cutoffs were not stringent enough. There also
# likely are other factors besides the ones we filtered for that are affecting the cell counts.

# After filtering for nGene per cell, you should still observe a small shoulder to the right of the
# main peak. What might this shoulder represent?
# Ans: This shoulder might represent doublets, which we did not address in our quality control
# process. Doublets are from when 2 cells were captured together, so they have the same barcode.
# This results in the gene counts from that barcode being higher than that of cells that were
# properly isolated, because you are actually counting genes detected from 2 cells instead of 1.

# When plotting the nGene against nUMI do you observe any data points in the bottom right quadrant
# of the plot? What can you say about these cells that have been removed?
# Ans: After filtering, there are no longer data points in the bottom right quadrant of the plot.
# These cells that were in the bottom right quadrant would have been cells with many transcripts
# but not very many genes (low novelty scores). As described in the tutorial, these cells could be
# cell types with low complexity, such as red blood cells, or they could just be contamination.
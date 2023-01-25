# Download files from NGDC and put them in this directory.
seuObject <- CreateSeuratObject(counts = data,  project = projectname,  min.cells = 3,  min.features = 200)
MySeurat.new <- function(RdsName,  patient,  pc,  res) {
  seu_object <- RdsName
  # filter
  mito.features <- grep(pattern = "^MT-",  x = rownames(x = seu_object),  value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seu_object,  slot = 'counts')[mito.features,  ]) / Matrix::colSums(x = GetAssayData(object = seu_object,  slot = 'counts'))
  seu_object[['percent.mito']] <- percent.mito
  
  # Save plots as pdf file
  pdf(file = paste(patient,  "_pc",  pc,  "_res",  res,  "_data_quality.pdf",  sep = ""),  width = 7,  height = 5.5)
  
  print(VlnPlot(object = seu_object,  features = c("nFeature_RNA",  "nCount_RNA",  "percent.mito"),  ncol = 3))
  print(FeatureScatter(object = seu_object,  feature1 = "nCount_RNA",  feature2 = "percent.mito"))
  print(FeatureScatter(object = seu_object,  feature1 = "nCount_RNA",  feature2 = "nFeature_RNA"))
  
  seu_object <- subset(x = seu_object,  subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.25)
  ## 2 Normalizing the data
  seu_object <- NormalizeData(object = seu_object,  normalization.method = "LogNormalize",  scale.factor = 1e4)
  ## 3 Detection of variable features across the single cell
  seu_object <- FindVariableFeatures(object = seu_object,  selection.method = 'mean.var.plot',  mean.cutoff = c(0.0125,  3),
                                     dispersion.cutoff = c(0.5,  Inf))
  length(x = VariableFeatures(object = seu_object))
  
  ### 4 Scaling the data and removing unwanted sources of variation
  # Remove technical noise and batch effects.
  seu_object <- ScaleData(object = seu_object, vars.to.regress = c("nCount_RNA",  "percent.mito"))
  
  ### 5 Perform linear dimensional reduction
  seu_object <- RunPCA(object = seu_object,  features = VariableFeatures(object = seu_object),  verbose = FALSE)
  
  ##### Determine statistically significant principal components
  print(ElbowPlot(object = seu_object, ndims = 50))
  dev.off()
  
  #### cluster cells
  seu_object <- FindNeighbors(object = seu_object,  dims = 1:pc)
  seu_object <- FindClusters(object = seu_object,  resolution = res)
  
  ### tSNE
  seu_object <- RunTSNE(object = seu_object,  dims = 1:pc)
  # note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
  
  pdf(file = paste(patient,  "_pc",  pc,  "_res",  res,  "_Dimplot.pdf",  sep = ""),  width = 7,  height = 5.5)
  print(DimPlot(object = seu_object,  reduction = 'tsne',  label = TRUE, raster = T))
  dev.off()
  # Save RDS files
  saveRDS(object = seu_object,  file = paste(patient,  "_pc",  pc,  "_res",  res,  ".rds",  sep = ""))
  
}
MySeurat.new(RdsName = seuObject,  patient = "All",  pc = 50,  res = 1)







#' Fig1
#' Fig1 B
#' tsne of all cells
#-----
library(Seurat)
library(ggsci)
library(ggplot2)

seuObject <- readRDS("All_patient_pc50_res1.rds")


### plot
# 1 set colors and orders
tissueorder <- c("Primary", "Portal-Venous", "Metastasis")
cellorder <- c("Epithelial", "Fibroblast", "Endothelial", "CTC", "B cell", "Myeloid", "NK", "T cell")

# tissue
library(scales)
mypal_tissue <- c("#62BFBF", "#D14A3D", "#62A1D9")
show_col(mypal_tissue)
names(mypal_tissue) <- tissueorder

# celltype
length(cellorder)
mypal_cell <- c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3", "#F781BE", "#33A9CF")
names(mypal_cell) <- cellorder
show_col(mypal_cell)

# patient
mypal_patient <- pal_npg()(6)

## 1.1 cell type
DimPlot(seuObject, label = T, group.by = "label2", repel = T, reduction = "tsne", cols = mypal_cell, 
        label.size = 4, label.color = mypal_cell, order = rev(cellorder), raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig1_All_dimplot_celltype.pdf", width = 6.5, height = 5)

## 1.2 tissue
DimPlot(seuObject, group.by = "group_tissue", reduction = "tsne", cols = mypal_tissue, raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig1_All_dimplot_tissue.pdf", width = 6.2, height = 5)


## 1.3 patient
DimPlot(seuObject, group.by = "group_patient", reduction = "tsne", cols = mypal_patient, raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig1_All_dimplot_patient.pdf", width = 5.7, height = 5)
#-----


#' Fig1
#' Fig1 C
#' tsne of tissue
#' (1) Blood
#' (2) Primary
#' (3) Metastasis
#-----
library(Seurat)
library(ggsci)
library(ggplot2)

### 1 Blood
seuObject <- readRDS(file = "All_Blood_pc40_res3.rds")

# 1 set colors and orders
cellorder <- c("T cell", "NK", "B cell", "Myeloid", "CTC")
library(scales)
# celltype
length(cellorder)
mypal_cell <- c("#33A9CF", "#F781BE", "#FF7F00", "#984EA3", "#E51A1D")
names(mypal_cell) <- cellorder
show_col(mypal_cell)

## 1.1 cell type
DimPlot(seuObject, label = T, group.by = "label2", repel = T, reduction = "tsne", cols = mypal_cell, 
        label.size = 4, label.color = mypal_cell, order = rev(cellorder), raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig1_All_blood_dimplot_celltype.pdf", width = 6.5, height = 5)

### 2 Primary
seuObject <- readRDS(file = "All_Primary_pc40_res1.rds")
table(seuObject$label2)

# 1 set colors and orders
cellorder <- c("Epithelial", "Fibroblast", "Endothelial", "B cell", "Myeloid", "NK", "T cell")

# colors
length(cellorder)
mypal_cell <- c("#377DB8", "#4DAF4A", "#A4572B", "#FF7F00", "#984EA3", "#F781BE", "#33A9CF")
names(mypal_cell) <- cellorder
show_col(mypal_cell)

## 1.1 cell type
DimPlot(seuObject, label = T, group.by = "label2", repel = T, reduction = "tsne", cols = mypal_cell, 
        label.size = 4, label.color = mypal_cell, order = rev(cellorder), raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig1_All_Primary_dimplot_celltype.pdf", width = 6, height = 5)


### 3 Metastasis
seuObject <- readRDS(file = "All_Metastasis_pc40_res1.rds")
table(seuObject$label2)

# 1 set colors and orders
cellorder <- c("Epithelial", "Fibroblast", "Endothelial", "B cell", "Myeloid", "NK", "T cell")

# colors
length(cellorder)
mypal_cell <- c("#377DB8", "#4DAF4A", "#A4572B", "#FF7F00", "#984EA3", "#F781BE", "#33A9CF")
names(mypal_cell) <- cellorder
show_col(mypal_cell)

## 1.1 cell type
DimPlot(seuObject, label = T, group.by = "label2", repel = T, reduction = "tsne", cols = mypal_cell, 
        label.size = 4, label.color = mypal_cell, order = rev(cellorder), raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig1_All_Metastasis_dimplot_celltype.pdf", width = 6, height = 5)
#-----




#' Fig1
#' Fig S3
#' cell type dotplot
#' (1) all cell
#' (2) HPV
#-----
library(Seurat)
library(ggsci)
library(ggplot2)


### 1 All cells
seuObject <- readRDS(file = "All_patient_pc50_res1.rds")

## 1.1 import data
genes <- c("EPCAM", "KRT8", "KRT18","KRT19",  # Epithelial
           "FAP", "COL1A1", "DCN", # Fibrobalst
           "VWF", "CDH5", "ENG", "PECAM1", # endothelial
           "PTPRC", "CD9", "TIMP1", "PPBP", "PF4", 
           "CD79A", "CD79B", "MS4A1",   # B cells
           "AIF1", "CD14", "LYZ", "FPR1",   # Myeloid
           "KLRF1", "KLRD1", "GNLY", # NK
           "CD3D", "CD3E", "CD3G") # T cells
cellorder <- c("Epithelial", "Fibroblast", "Endothelial", "CTC", "B cell", "Myeloid", "NK", "T cell")
# groups of marker gene
geneorder <- c("Epithelial", "Fibroblast", "Endothelial", "CTC", "B cell", "Myeloid", "NK", "T cell")
groups <- rep(geneorder, c(4, 3, 4, 5, 3, 4, 3, 3))

## plot
cellDotplot <- function(seuObject, genes, cellorder, geneorder, groups) {
  library(tidyverse)
  library(ggplot2)
  count <- as.matrix(GetAssayData(object = seuObject, slot = "data"))[genes, , drop = FALSE]
  count <- t(count)
  count[1:5, 1:5]
  
  # 1.1 mean expression of genes in every clusters
  group <- as.character(seuObject$label)
  gene.mean <- rowsum(count, group = group) / as.vector(table(group))
  #gene.mean <- scale(gene.mean)
  gene.mean[gene.mean > 1.5] <- 1.5
  gene.mean <- as.data.frame(t(gene.mean))
  gene.mean$gene <- rownames(gene.mean)
  gene.mean <- pivot_longer(gene.mean, !gene, names_to = "label", values_to = "z_score")
  gene.mean$group <- rep(groups, each = length(cellorder))
  
  # 1.2 percent of expressed cells
  percent <- count
  percent[percent > 0] <- 1
  gene.percent <- rowsum(percent, group = group) / as.vector(table(group))
  gene.percent <- as.data.frame(t(gene.percent))
  
  gene.percent$gene <- rownames(gene.percent)
  gene.percent <- pivot_longer(gene.percent, !gene, names_to = "label", values_to = "percent")
  
  # 1.3 combine
  df.plot <- cbind(gene.mean, percent = gene.percent$percent)
  
  ## 2 plot
  p <- df.plot %>% 
    ggplot(aes(factor(gene, levels = genes), 
               factor(label, levels = rev(cellorder)))) + 
    geom_point(aes(size = percent, colour = z_score)) + 
    theme(strip.text.x = element_blank(),  
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 13), 
          legend.title = element_text(size = 13), 
          legend.text = element_text(size = 13), 
          axis.text.y = element_text(color = "black"), 
          axis.text.x = element_text(color = "black", angle = -90, hjust = 0), 
          panel.background = element_rect(colour = "black", fill = "white"), 
          panel.grid = element_line(colour = "grey", linetype = "dashed"), 
          panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)) + 
    facet_grid(. ~ factor(group, levels = geneorder), scales = "free", space = "free") + 
    scale_colour_distiller(palette = "RdYlBu") +
    labs(x = "", y = "")
  return(p)
}

cellDotplot(seuObject = seuObject, genes = genes, cellorder = cellorder, geneorder = geneorder, groups = groups)
ggsave(filename = "Fig1_All_celltype_dotplot_expression_all_cell.pdf", width = 8.5, height = 5)


#### 2 P12345 blood
seu.blood <- readRDS(file = "All_Blood_pc40_res3.rds")

genes <- c("CD3D", "CD3E", "CD3G", # T cells
           "KLRF1", "KLRD1", "FGFBP2", "NKG7", # NK
           "CD79A", "CD79B", "MS4A1",   # B cells
           "AIF1", "CD14", "LYZ", "FPR1",  # meyloid
           "PTPRC", "CD9", "KRT8", "MMP1", "TIMP1") 
cellorder <- c("T cell", "NK", "B cell", "Myeloid", "CTC")
# groups of marker gene
geneorder <- c("T_cell", "NK", "B_cell", "Myeloid", "CTC")
groups <- rep(geneorder, c(3, 4, 3, 4, 5))

## plot
cellDotplot(seuObject = seu.blood, genes = genes, cellorder = cellorder, geneorder = geneorder, groups = groups)
ggsave(filename = "Fig1_All_Blood_cell_type_dotplot.pdf", width = 6, height = 3)
#-----




#' Fig1
#' Fig S3
#' cell type
#' featureplot
#-----
library(Seurat)
library(ggsci)
library(ggplot2)

seuObject <- readRDS(file = "../P12345_pc50_res1_plot.rds")
genes <- c("EPCAM", "KRT18", # Epithelial
           "FAP", "COL1A1", # Fibrobalst
           "VWF", "PECAM1", # endothelial
           "PTPRC", "PPBP",  
           "CD79A", "CD19",   # B cells
           "AIF1", "CD14", "LYZ", "FPR1",   # Myeloid
           "KLRF1", "KLRD1", # NK
           "CD3D", "CD3G")
FeaturePlot(seuObject, features = genes, reduction = "tsne", raster = T, ncol = 6)
ggsave(filename = "Fig1_celltype_featureplot_genes.pdf", width = 24, height = 10.5)
#-----




#' Fig1
#' Fig S3
#' Primary+Blood+Metastasis
#' barplot
#-----
library(Seurat)

setwd("/home/SJE/project/pancreatic_scRNA/analysis/Figures/2022_P6/")
p12345 <- readRDS(file = "20220323_all_patient_pc50_res1.rds")
table(p12345$label2)
label <- p12345$label2

# tissue
tissue <- as.character(p12345$group_tissue)
table(tissue)
tissue[tissue == "Portal-Venous"] <- "Blood"
p12345$group_tissue <- tissue

seuObject <- p12345
tissueorder <- c("Primary", "Blood", "Metastasis")
table(seuObject$label2)

# colors
library(scales)
mypal_tissue <- c("#62BFBF", "#D14A3D", "#62A1D9")
show_col(mypal_tissue)
names(mypal_tissue) <- tissueorder

cellorder <- c("Epithelial", "Fibroblast", "Endothelial", "B cell", "Myeloid", "NK", "T cell", "CTC")
length(cellorder)
mypal_cell <- c("#377DB8", "#4DAF4A", "#A4572B", "#FF7F00", "#984EA3", "#F781BE", "#33A9CF", "#E51A1D")
names(mypal_cell) <- cellorder
show_col(mypal_cell)

## 1.1 Group by cell type
df <- data.frame("patient" = seuObject$group_patient, 
                 "tissue" = factor(seuObject$group_tissue, levels = tissueorder), 
                 "celltype" =  factor(seuObject$label2, levels = cellorder))
df.group <- group_by(df, tissue, patient, celltype)
ratio <- summarise(df.group, num = n())
cellgroup <- as.character(ratio$celltype)
cellgroup[cellgroup %in% c("Epithelial", "Fibroblast", "Endothelial")] <- "Non immune cells"
cellgroup[cellgroup %in% c("B cell", "Myeloid", "NK", "T cell", "CTC")] <- "Immune cells"
table(cellgroup)
ratio$cellgroup <- factor(cellgroup, levels = c("Non immune cells", "Immune cells"))

ratio <- group_by(ratio, tissue, patient, cellgroup) %>% mutate(ratio.group = num / sum(num))

library(ggplot2)
library(ggpubr)
df.plot <- ratio
df.plot$name <- factor(paste(df.plot$patient, df.plot$tissue, sep = " "), 
                       levels = unique(paste(df.plot$patient, df.plot$tissue, sep = " ")))
head(df.plot)

ggplot(df.plot, aes(x=cellgroup, y=ratio.group, fill=celltype, group=cellgroup))+
  geom_bar(stat="identity",position="stack") +
  theme_classic() +
  facet_grid(.~name) +
  ylab("Ratio of cell type") +
  scale_fill_manual(values = mypal_cell) +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), 
        axis.title = element_text(color = "black", size = 10), 
        axis.text.y = element_text(color = "black", size = 10), 
        axis.text.x = element_blank(),
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 10, angle = 0))
ggsave(filename = "Fig1_barplot_fraction_of_cells.pdf", width = 5.6, height = 3.5)
#-----




#' Fig1
#' Fig S4
#' Epithelial cells + CTCs
#-----
library(Seurat)
MySeurat.new <- function(RdsName,  patient,  pc,  res) {
  seu_object <- RdsName
  # filter
  mito.features <- grep(pattern = "^MT-",  x = rownames(x = seu_object),  value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seu_object,  slot = 'counts')[mito.features,  ]) / Matrix::colSums(x = GetAssayData(object = seu_object,  slot = 'counts'))
  seu_object[['percent.mito']] <- percent.mito
  
  # Save plots as pdf file
  pdf(file = paste(patient,  "_pc",  pc,  "_res",  res,  "_data_quality.pdf",  sep = ""),  width = 7,  height = 5.5)
  
  print(VlnPlot(object = seu_object,  features = c("nFeature_RNA",  "nCount_RNA",  "percent.mito"),  ncol = 3))
  print(FeatureScatter(object = seu_object,  feature1 = "nCount_RNA",  feature2 = "percent.mito"))
  print(FeatureScatter(object = seu_object,  feature1 = "nCount_RNA",  feature2 = "nFeature_RNA"))
  
  seu_object <- subset(x = seu_object,  subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.25)
  ## 2 Normalizing the data
  seu_object <- NormalizeData(object = seu_object,  normalization.method = "LogNormalize",  scale.factor = 1e4)
  ## 3 Detection of variable features across the single cell
  seu_object <- FindVariableFeatures(object = seu_object,  selection.method = 'mean.var.plot',  mean.cutoff = c(0.0125,  3),
                                     dispersion.cutoff = c(0.5,  Inf))
  length(x = VariableFeatures(object = seu_object))
  
  ### 4 Scaling the data and removing unwanted sources of variation
  # Remove technical noise and batch effects.
  seu_object <- ScaleData(object = seu_object, vars.to.regress = c("nCount_RNA",  "percent.mito"))
  
  ### 5 Perform linear dimensional reduction
  seu_object <- RunPCA(object = seu_object,  features = VariableFeatures(object = seu_object),  verbose = FALSE)
  
  ##### Determine statistically significant principal components
  print(ElbowPlot(object = seu_object, ndims = 50))
  dev.off()
  
  #### cluster cells
  seu_object <- FindNeighbors(object = seu_object,  dims = 1:pc)
  seu_object <- FindClusters(object = seu_object,  resolution = res)
  
  ### tSNE
  seu_object <- RunTSNE(object = seu_object,  dims = 1:pc)
  # note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
  
  pdf(file = paste(patient,  "_pc",  pc,  "_res",  res,  "_Dimplot.pdf",  sep = ""),  width = 7,  height = 5.5)
  print(DimPlot(object = seu_object,  reduction = 'tsne',  label = TRUE, raster = T))
  dev.off()
  # Save RDS files
  saveRDS(object = seu_object,  file = paste(patient,  "_pc",  pc,  "_res",  res,  ".rds",  sep = ""))
  
}

pall <- readRDS(file = "All_patient_pc50_res1.rds")
table(pall$label2)

## 1 Epi + CTC
seuObject <- subset(pall, subset=label2 %in% c("Epithelial", "CTC"))
MySeurat.new(RdsName = seuObject,  patient = "All_epi_ctc",  pc = 40,  res = 1)
#-----



#' Fig1
#' Fig 1D
#' DEgenes
#' CTC/P+M tumor cells
#' top30
#' heatmap
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

seuObject <- readRDS(file = "All_all_epi_ctc_pc40_res1.rds")
table(seuObject$label)

# 去除normal epithelial cells
seuObject <- subset(seuObject, subset=label != "Normal epithelial")
Idents(seuObject) <- "label"

### 2 DE genes
degenes <-FindMarkers(seuObject, ident.1 = "CTC", "Tumor cells", logfc.threshold = 0.01, only.pos = F)
degenes <- degenes[order(degenes$avg_log2FC, decreasing = T), ]


seu.sub <- seuObject

### 1 Plot : heatmap
group <- as.character(seu.sub$group_tissue)
table(group)
group[group=="Portal-Venous"] <- "CTC"

order1 <- c("Primary", "CTC", "Metastasis")
orders <- c("CTC", "Epithelial")
# colors
cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

mypal.tissue <- c("#62BFBF", "#D14A3D", "#62A1D9")
names(mypal.tissue) <- order1

# annotation
ha.t <- HeatmapAnnotation(Group = factor(group, levels = order1),
                          col = list(Group = mypal.tissue))

# genes
top = 30
degenes$gene <- rownames(degenes)
topn <- degenes[c(1:top, (nrow(degenes)-(top-1)):nrow(degenes)), ]
topn$cluster <- rep(orders, c(top, top))


# data to plot
ht_opt$message = FALSE
count <- as.matrix(GetAssayData(object = seu.sub, slot = "data"))[rownames(topn), , drop = FALSE]
exMatrix <- t(scale(t(count)))

# plot
width = 9
height <- nrow(topn) * 0.2 + 2
pdf("Fig1_D_CTC_vs_epithelial_Degenes_top30.pdf", height = height, width = width)
Heatmap(exMatrix,
        show_column_names = F,
        show_row_names = T,
        row_title = NULL,
        column_split = factor(group, levels = order1),
        row_split = factor(topn$cluster, levels = orders),
        cluster_columns = T,
        cluster_rows = F,
        show_column_dend = F,
        cluster_column_slices = F,
        top_annotation = ha.t,
        col = cols,
        name = "Z score",
        column_title_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 15),
        row_names_max_width = max_text_width(rownames(exMatrix),
                                             gp = gpar(fontsize = 15)))
dev.off()
#-----




#' Fig1
#' Fig 1E
#' enrichment analysis
#' CTC
#' Heatmap
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggsci)
library(ComplexHeatmap)
library(circlize)


seuObject <- readRDS(file = "All_all_epi_ctc_pc40_res1.rds")
table(seuObject$label)

tissue <- "Pall_tumorcell_ctc"

## 1 subset
cells <- c("Tumor cells", "CTC")
seu.sub <- subset(seuObject, subset=label %in% cells)
# Output data for Single Cell Signature Score anslyais
library(data.table)
exMatrix <- as.matrix(GetAssayData(seu.sub, slot = "data"))
exMatrix <- exMatrix[rowSums(exMatrix) > 10, ]
exMatrix <- exMatrix[rowSums(exMatrix > 0) > 10, ]

exMatrix[1:5, 1:5]
exMatrix <- round(exMatrix, digits = 2)
exMatrix <- t(exMatrix)
# set as data.frame for row names
exMatrix <- as.data.frame(exMatrix)
setwd("~/software/SingleCellSignatureExplorer/SingleCellSignatureScorer/data")
fwrite(exMatrix, file = "Pall_tumorcell_ctc_data.tsv", sep = "\t", row.names = T, col.names = T)
## Process Single Cell Signature Score anslyais in Shell.


## 2 input scss score
score1 <- fread(file = paste0("../../scss_output/pall_tumor_ctc/H_", tissue,"_data.tsv"), header = T, data.table = F)
score2 <- fread(file = paste0("../../scss_output/pall_tumor_ctc/C2_CP_KEGG_", tissue,"_data.tsv"), header = T, data.table = F)
score3 <- fread(file = paste0("../../scss_output/pall_tumor_ctc/C5_BP_", tissue,"_data.tsv"), header = T, data.table = F)
score4 <- fread(file = paste0("../../scss_output/pall_tumor_ctc/C6_", tissue,"_data.tsv"), header = T, data.table = F)
score5 <- fread(file = paste0("../../scss_output/pall_tumor_ctc/C2_CP_REACTOME_", tissue,"_data.tsv"), header = T, data.table = F)
gsvascore <- cbind(score1, score2, score3, score4, score5)
gsvascore[1:5, 1:5]

## 3 geneset to plot
genesets <- read.csv(file = "../../data/CTC_genesets_heatmap", header = T)
table(genesets$Genesets %in% colnames(gsvascore))
genesets <- genesets[genesets$Genesets %in% colnames(gsvascore), ]

df.plot <- gsvascore[, genesets$Genesets]
rownames(df.plot) <- gsvascore$id
df.plot <- t(df.plot)
df.plot[1:5, 1:5]
df.plot <- df.plot[, colnames(seu.sub)]

df.plot <- t(scale(t(df.plot)))

## 4 plot
group <- as.character(seu.sub$group_tissue)
group[group == "Portal-Venous"] <- "CTC"
table(group)
orders <- c("Primary", "CTC", "Metastasis")
group <- factor(group, levels = orders)

# colors
cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

mypal.tissue <- c("#62BFBF", "#D14A3D", "#62A1D9")
names(mypal.tissue) <- orders

# annotation
ha.t <- HeatmapAnnotation(Group = group,
                          col = list(Group = mypal.tissue))
x.length <- unlist(lapply(strsplit(as.character(rownames(df.plot)), split = ""), length))
width = (max(x.length) * 0.1) + 6
height = (0.1 * nrow(df.plot) + 1)

genesets$Group <- factor(genesets$group, levels = unique(genesets$group))
group.row <- as.data.frame(table(genesets$Group))
col <- c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3", "#F781BE", "#33A9CF")
row.color <- rep(col, group.row$Freq)

## 1 plot
pdf("Fig1_E_P123456_CTC_pathays_heatmap.pdf", height = height, width = width)
Heatmap(df.plot, 
        cluster_rows = F, 
        cluster_row_slices = genesets$group,
        show_row_dend = F, 
        cluster_columns = F, 
        row_title = NULL, 
        column_title = NULL,
        column_split = group, 
        row_split = factor(genesets$Group, levels = unique(genesets$Group)), 
        top_annotation = ha.t, 
        col = cols, 
        name = "Z score", 
        show_column_names = F,
        column_title_gp = gpar(fontsize = 5), 
        row_names_gp = gpar(fontsize = 6, col = row.color), 
        row_names_max_width = max_text_width(rownames(df.plot), 
                                             gp = gpar(fontsize = 6)))
dev.off()

#-----





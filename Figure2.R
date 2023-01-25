# Download files from NGDC and put them in this directory.
# Immune cells only



#' Figure 2
#' Figure 2A-2B
#' All immune cells
#' Dimplot
#-----
library(Seurat)
library(ggplot2)

seu.immu <- readRDS("All_immune_pc40_res1.rds")

### 1 Dimplot
tissueorder <- c("Portal-Venous", "Primary", "Metastasis")
cellorder <- c("B cell", "Naive T", "CD8 EFF", "CD8 Ex", "Memory T", "Treg", "NKT", "NK", 
               "cDC", "pDC", "M1", "M2", "Monocyte", "Neutrophil", "Mast")

# tissue
library(scales)
mypal_tissue <- c("#D14A3D", "#62BFBF", "#62A1D9")
show_col(mypal_tissue)
names(mypal_tissue) <- tissueorder

# celltype
length(cellorder)
mypal_cell <- c(pal_npg()(7), 
                c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3",
                  "#F781BE", "#33A9CF"))
show_col(mypal_cell)
names(mypal_cell) <- cellorder
show_col(mypal_cell)


## 1.1 cell type
DimPlot(seu.immu, label = T, group.by = "label", repel = T, reduction = "tsne", cols = mypal_cell, 
        label.size = 4, order = rev(cellorder), raster = T, label.color = mypal_cell) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig2_All_immune_dimplot_celltype.pdf", width = 6.3, height = 5)

DimPlot(seu.immu, label = F, group.by = "group_tissue", repel = T, reduction = "tsne", cols = mypal_tissue, 
        raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig3_P12345Fig2_All_immune_dimplot_tissue.pdf", width = 6.3, height = 5)
#-----






#' Figure 2
#' Figure 2C
#' Plot of mark genes (immune cells)
#' heatmap
#-----
library(Seurat)
library(ggsci)
library(ggplot2)


### 1
seuObject <- readRDS(file = "All_immune_pc40_res1.rds")
table(seuObject$label)
Idents(seuObject) <- seuObject$label

## 1.1 import data
genes <- c("KLRB1","KLRF1", "KLRD1", "NCAM1","CD3D", "CD3E", "CD3G",  # NK
           "CD8A","LAG3", "TIGIT", "CTLA4", "HAVCR2", "PDCD1", # CD8 Ex
           "GZMA", "GZMB","IFNG","GZMK", # CD8 EFF
           "FOXP3", "IL2RA", "IKZF2",   # Treg
           "CD44", "IL7R", "LTB", #Memory T
           "CCR7", "TCF7", "LEF1","SELL", # Naive T
           "CD79A", "CD79B", "MS4A1", # B cell
           "FCGR3B", "FPR1", # Neutrophil
           "CD14", "S100A12", "FCGR3A", # Monocyte
           "CD68", "ITGAX", "ITGAM", "CD86","IL1B",  # M1
           "CD163", "MRC1", "MSR1", # M2
           "CD1C", "FCER1A", "CLEC10A", # cDC
           "LILRA4","CCDC50","IL3RA",  # pDC
           "MS4A2", "TPSAB1", "KIT")  # Mast 

### 2 Hetamap
count <- as.matrix(GetAssayData(object = seuObject, slot = "data"))[genes, , drop = FALSE]
count <- as.data.frame(t(count))
group <- as.character(Idents(seuObject))
gene.mean <- rowsum(count, group = group) / as.vector(table(group))
gene.mean <- scale(gene.mean)
gene.mean <- t(gene.mean)[, cellorder]

library(ComplexHeatmap)
library(circlize)

cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

col.split <- rep(c("NK/NKT", "T", "B", "Myeloid"), c(2, 5, 1, 7))
ha.t <- HeatmapAnnotation(Group = factor(col.split, levels = unique(col.split)))
mypal_cell <- c(pal_npg()(6), 
                c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3",
                  "#F781BE", "#33A9CF"))
show_col(mypal_cell)
order1 <- c("B cell", "Naive T", "CD8 EFF", "CD8 Ex", "Memory T", "Treg", "NK", 
            "cDC", "pDC", "M1", "M2", "Monocyte", "Neutrophil", "Mast")
names(mypal_cell) <- order1
mypal_cell <- mypal_cell[geneorder]
show_col(mypal_cell)
ha.l <- rowAnnotation("Cell type" = factor(groups, levels = geneorder), 
                      col=list("Cell type"=mypal_cell))
pdf(file = "Fig2_D_All_immune_celltype_Heatmap.pdf", width = 9.2, height = 9.5)
Heatmap(gene.mean, 
        cluster_rows = F, cluster_columns = F,
        col = cols, 
        row_split = factor(groups, levels = geneorder), 
        column_split = factor(cellorder, levels = unique(cellorder)),
        name = "Z score", 
        column_title = NULL,
        column_names_side = "top", 
        column_names_rot = 45, 
        row_title_rot = 90,
        top_annotation = ha.t, 
        left_annotation = ha.l,
        layer_fun = function(j, i, x, y, width, height, fill) {
          v = pindex(gene.mean, i, j)
          if(sum(v > 0)/length(v) > 0.75) {
            grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
          }
        }
)
dev.off()
#-----






#' Figure 2
#' Figure S7
#' Plot of mark genes (immune cells)
#' featureplot
#-----
library(Seurat)
library(ggsci)
library(ggplot2)


### 1 import data
seuObject <- readRDS(file = "All_immune_pc40_res1.rds")
### 5 featrueplot
genes <- c("KLRF1", "KLRD1", "GNLY","CD3D", "CD3E","CD3G",  # NK
           "CD8A","LAG3", "TIGIT", "CTLA4", "HAVCR2", "PDCD1", # CD8 Ex
           "GZMA", "IFNG","GZMK", # CD8 EFF
           "FOXP3", "IL2RA", "IKZF2",   # Treg
           "CD44", "IL7R", "LTB", #Memory T
           "CCR7", "TCF7", "LEF1",# Naive T
           "CD79A", "CD79B", "MS4A1", # B cell
           "FCGR3B", "FPR1","MNDA", # Neutrophil
           "CD14", "S100A12", "FCGR3A", # Monocyte
           "CD68", "ITGAX", "ITGAM",  # M1
           "CD163", "MRC1", "MSR1", # M2
           "CD1C", "FCER1A", "CLEC10A", # cDC
           "LILRA4","CCDC50","IL3RA",  # pDC
           "MS4A2", "TPSAB1", "KIT")  # Mast 
FeaturePlot(seuObject, features = genes, reduction = "tsne", ncol = 6, raster = T)+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(filename = "Figure3_p12345_immune_celltype_featureplot_all_gene.pdf", width = 24, height = 28)
#-----





#' Figure 2
#' Figure 2D-2G
#' CellPhoneDB
#-----
library(Seurat)
library(data.table)

### input
p123456 <- readRDS(file = "All_patient_pc50_res1.rds")
seu.immu <- readRDS("All_immune_pc40_res1.rds")
table(seu.immu$label)

seu.m <- subset(p123456, subset=label %in% c("Tumor cells", "CTC"))
p123456 <- merge(seu.immu, seu.m)

#### 1 Primary: 
seu.sub <- subset(p123456, subset = group_tissue %in% c("Primary"))

# cell type
group <- seu.sub$label
table(group)
group[group=="Tumor cells"] <- "Primary"

# combine
df <- as.data.frame(as.matrix(GetAssayData(seu.sub, slot = "data")))
df <- df[rowSums(df) > 10, ]
Gene <- rownames(df)
df <- cbind(Gene, df)
df[1:3, 1:3]

celltype <- data.frame("Cell" = colnames(seu.sub), "Cell_type" = group)
table(celltype$Cell_type)
dir.create("Primary")
fwrite(df, file = "Primary/cellphonedb_count.txt", row.names = F, col.names = T, sep = "\t", quote = F)
fwrite(celltype, file = "Primary/cellphonedb_celltype.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#### 2 Metastasis: 
seu.sub <- subset(p123456, subset = group_tissue %in% c("Metastasis"))

# cell type
group <- seu.sub$label
table(group)
group[group=="Tumor cells"] <- "Metastasis"

# combine
df <- as.data.frame(as.matrix(GetAssayData(seu.sub, slot = "data")))
df <- df[rowSums(df) > 10, ]
Gene <- rownames(df)
df <- cbind(Gene, df)
df[1:5, 1:5]

celltype <- data.frame("Cell" = colnames(seu.sub), "Cell_type" = group)
table(celltype$Cell_type)
dir.create("Metastasis")
fwrite(df, file = "Metastasis/cellphonedb_count.txt", row.names = F, col.names = T, sep = "\t", quote = F)
fwrite(celltype, file = "Metastasis/cellphonedb_celltype.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#### 3 Blood: 
seu.sub <- subset(p123456, subset = group_tissue %in% c("Portal-Venous"))

# cell type
group <- seu.sub$label
table(group)

# combine
df <- as.data.frame(as.matrix(GetAssayData(seu.sub, slot = "data")))
df <- df[rowSums(df) > 10, ]
Gene <- rownames(df)
df <- cbind(Gene, df)
df[1:5, 1:5]

celltype <- data.frame("Cell" = colnames(seu.sub), "Cell_type" = group)
dir.create("Blood")
fwrite(df, file = "Blood/cellphonedb_count.txt", row.names = F, col.names = T, sep = "\t", quote = F)
fwrite(celltype, file = "Blood/cellphonedb_celltype.txt", row.names = F, col.names = T, sep = "\t", quote = F)


saveRDS(p123456, file = "20220611_Pall_tumor_ctc_immune.rds")

### 4 CellPhoneDB analysis in shell
# cd ~/project/pancreatic_scRNA/analysis/Figures/cellphonedb/Primary/
  
# cellphonedb method statistical_analysis cellphonedb_celltype.txt \
# cellphonedb_count.txt --iterations=100 \
# --threads=12 --counts-data=gene_name --output-path=outs
#-----
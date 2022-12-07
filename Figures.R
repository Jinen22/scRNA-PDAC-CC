#' 2021/06/27
#' Fig2
#' Fig 2.1
#' tsne of non immune cells
#-----
library(Seurat)
library(ggsci)
library(ggplot2)

setwd("/home/SJE/project/pancreatic_scRNA/analysis/Figures/Figure2/")
p12345.non <- readRDS("../P12345_non_immune_pc40_res0.5.rds")
table(p12345.non$group_tissue)
table(p12345.non$label2)
table(p12345.non$label)

# 1 set colors and orders
tissueorder <- c("Portal-Venous", "Primary", "Metastasis")
cellorder <- c("Epithelial", "Fibroblast", "Endothelial", "CTC")

# tissue
library(scales)
mypal_tissue <- c("#D14A3D", "#62BFBF", "#62A1D9")
show_col(mypal_tissue)
names(mypal_tissue) <- tissueorder

# celltype
length(cellorder)
mypal_cell <- c("#377DB8", "#4DAF4A", "#F8DB35", "#E51A1D")
names(mypal_cell) <- cellorder
show_col(mypal_cell)


# patient
mypal_patient <- pal_npg()(5)

## 1.1 cell type
DimPlot(p12345.non, label = T, group.by = "label2", repel = T, reduction = "tsne", cols = mypal_cell, 
        label.size = 4, label.color = mypal_cell, order = rev(cellorder), raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig2_P12345_non_immune_dimplot_celltype.pdf", width = 6.2, height = 5)

## 1.2 tissue
DimPlot(p12345.non, group.by = "group_tissue", reduction = "tsne", cols = mypal_tissue, raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig2_P12345_non_immune_dimplot_tissue.pdf", width = 6.2, height = 5)


## 1.3 patient
DimPlot(p12345.non, group.by = "group_patient", reduction = "tsne", cols = mypal_patient, raster = T) +
  theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "black", size = 16), 
        axis.title = element_text(colour = "black", size = 16))
ggsave(filename = "Fig2_P12345_non_immune_dimplot_patient.pdf", width = 5.7, height = 5)



FeaturePlot(p12345.non, features = c("PTPRC", "EPCAM", "CDH1", "PDX1"))

DimPlot(p12345, group.by = "group_patient", reduction = "tsne")
DimPlot(p12345, label = T, group.by = "label")
FeaturePlot(seu.sub, features = c("PTPRC", "EPCAM", "CDH1", "PDX1"))
seu.sub <- subset(p12345, subset = label == "Islets")
DotPlot(seu.sub, features = "PDX1")
DimPlot(seu.sub, group.by = "group_tissue")
table(as.character(seu.sub$seurat_clusters))
seu.sub <- subset(p12345, idents = 43)
table(seu.sub$group_tissue)
#-----





#' 2021/07/01
#' Fig2
#' 相关性图分析
#-----
library(Seurat)
library(ggsci)
library(ggplot2)

setwd("/home/SJE/project/pancreatic_scRNA/analysis/Figures/Figure2/")
p12345 <- readRDS(file = "../P12345_pc50_res1_plot.rds")

seu.sub <- subset(p12345, subset = label2 %in% c("Epithelial", "CTC"))
df <- as.matrix(GetAssayData(seu.sub, slot = "data"))
df <- df[rowSums(df) > 100, ]

### 1 合并在一起
group <- as.character(seu.sub$group_tissue)
table(group)
group[group == "Portal-Venous"] <- "CTC"
gene.mean <- rowsum(t(df), group = group) / as.vector(table(group))
gene.mean <- t(gene.mean)

library(corrplot)
matrix <- cor(gene.mean,method="pearson")
orders <- c("CTC", "Priamry", "Metastasis")

cols = colorRamp2(seq(from = 0, to = 1, length.out = 5),
                  c("#FDDBC7", "#F4A582", "#D6604D","#B2182B", "#67001F"))
mypal_tissue <- c("#D14A3D", "#62BFBF", "#62A1D9")
show_col(mypal_tissue)
names(mypal_tissue) <- orders
ha.t <- HeatmapAnnotation(Tissue = factor(orders, levels = orders),col = list(Tissue=mypal_tissue))
ha.l <- rowAnnotation(Tissue = factor(orders, levels = orders), col=list(Tissue=mypal_tissue))

pdf("20220105_FigS5_corrplot_CTC_epi.pdf", width = 5, height = 4)
Heatmap(matrix, cluster_rows = F,  cluster_columns = F, show_row_names = F, show_column_names = F, 
        col = cols, 
        row_split = orders, 
        top_annotation = ha.t, 
        left_annotation = ha.l,
        column_split = orders,
        name = "Correlation score")
dev.off()


# 1 set colors and orders

### 2 所有细胞
tissueorder <- c("Primary", "Portal-Venous", "Metastasis")

# matrix <- cor(df,method="pearson")
matrix <- fread("CTC_M_P_coorplot.csv", data.table = F)
rownames(matrix) <- matrix$V1
matrix <- matrix[, -1]
matrix[1:3, 1:3]
group <- seu.sub$group_tissue
group <- factor(as.character(group), levels = tissueorder)
orders <- order(group)
df.use <- matrix[orders, orders]
group <- group[orders]

df.use[1:3, 1:3]

library(ComplexHeatmap)
png("20220124_Fig2_CTC_P_M_corplot_all_cells.png", width = 950, height = 800)
Heatmap(df.use, cluster_rows = F,  cluster_columns = F, 
        show_row_names = F, show_column_names = F, 
        col = cols, 
        row_split = group, 
        column_split = group, 
        name = "Score")
dev.off()
#-----


#' 2021/07/16
#' DEgenes
#' P12345 epithelial : Primary / Metastasis
#' top50
#' (1) heatmap
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

dirs <- "~/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)

seuObject <- readRDS(file = "../P12345_pc50_res1_plot.rds")
table(seuObject$label)
Idents(seuObject) <- "label"

### 2 DE genes
degenes <- read.csv(file = "../../DEgenes/P_M/P12345_eip_Primary_vs_Metastasis_degenes_2th.csv", header = T, row.names = 1)
degenes <- degenes[order(degenes$avg_log2FC, decreasing = T), ]

library(ggsci)
library(ComplexHeatmap)
library(circlize)

cells <- c("Epithelial")
seu.sub <- subset(seuObject, subset = label %in% cells)

### 1 Plot : heatmap
group <- as.character(seu.sub$group_tissue)

# plot : heatmap
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
table(seu.sub$group_tissue)

order1 <- c("Primary", "Metastasis")
orders <- c("Primary", "Metastasis")
# colors
cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

mypal.tissue <- c("#62BFBF", "#62A1D9")
names(mypal.tissue) <- order1

# annotation
ha.t <- HeatmapAnnotation(Group = factor(group, levels = order1),
                          col = list(Group = mypal.tissue))

# genes
degenes <- degenes[order(degenes$avg_log2FC, decreasing = T), ]
top = 50
topn <- degenes[c(1:50, (nrow(degenes)-49):nrow(degenes)), ]

# data to plot
ht_opt$message = FALSE
count <- as.matrix(GetAssayData(object = seu.sub, slot = "data"))[rownames(topn), , drop = FALSE]
exMatrix <- t(scale(t(count)))

# plot
width = 10
height <- nrow(topn) * 0.1 + 2
pdf("Fig2_P12345_Primary_vs_metastasis.pdf", height = height, width = width)
Heatmap(exMatrix,
        show_column_names = F,
        show_row_names = T,
        row_title = NULL,
        column_split = factor(group, levels = order1),
        row_split = factor(rep(orders, c(50, 50)), levels = orders),
        cluster_columns = T,
        cluster_rows = F,
        show_column_dend = F,
        cluster_column_slices = F,
        top_annotation = ha.t,
        col = cols,
        name = "Z score",
        column_title_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        row_names_max_width = max_text_width(rownames(exMatrix),
                                             gp = gpar(fontsize = 10)))
dev.off()


#-----




#' 2021/07/16
#' DEgenes
#' P12345 epithelial : Primary / Metastasis
#' top50
#' (1) heatmap
#' (2) volcano plot
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

dirs <- "~/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)

seuObject <- readRDS(file = "../P12345_pc50_res1_plot.rds")
table(seuObject$label)
Idents(seuObject) <- "label"

### 2 DE genes
degenes <- read.csv(file = "../../DEgenes/P_M/P12345_eip_Primary_vs_Metastasis_degenes_2th.csv", header = T, row.names = 1)
degenes <- degenes[order(degenes$avg_log2FC, decreasing = T), ]

library(ggsci)
library(ComplexHeatmap)
library(circlize)

cells <- c("Epithelial")
seu.sub <- subset(seuObject, subset = label %in% cells)

### 1 Plot : heatmap
group <- as.character(seu.sub$group_tissue)

# plot : heatmap
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
table(seu.sub$group_tissue)

order1 <- c("Primary", "Metastasis")
orders <- c("Primary", "Metastasis")
# colors
cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

mypal.tissue <- c("#62BFBF", "#62A1D9")
names(mypal.tissue) <- order1

# annotation
ha.t <- HeatmapAnnotation(Group = factor(group, levels = order1),
                          col = list(Group = mypal.tissue))

# genes
degenes <- degenes[order(degenes$avg_log2FC, decreasing = T), ]
top = 50
topn <- degenes[c(1:50, (nrow(degenes)-49):nrow(degenes)), ]

# data to plot
ht_opt$message = FALSE
count <- as.matrix(GetAssayData(object = seu.sub, slot = "data"))[rownames(topn), , drop = FALSE]
exMatrix <- t(scale(t(count)))

# plot
width = 10
height <- nrow(topn) * 0.1 + 2
pdf("Fig2_P12345_Primary_vs_metastasis.pdf", height = height, width = width)
Heatmap(exMatrix,
        show_column_names = F,
        show_row_names = T,
        row_title = NULL,
        column_split = factor(group, levels = order1),
        row_split = factor(rep(orders, c(50, 50)), levels = orders),
        cluster_columns = T,
        cluster_rows = F,
        show_column_dend = F,
        cluster_column_slices = F,
        top_annotation = ha.t,
        col = cols,
        name = "Z score",
        column_title_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        row_names_max_width = max_text_width(rownames(exMatrix),
                                             gp = gpar(fontsize = 10)))
dev.off()

### 2 volcano plot
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)

### 1 input
#只保留差异倍数和P值这两列
head(degenes)
tail(degenes)
degenes <- degenes[, c(2, 5)]

# reneme names
colnames(degenes) <- c("avg_logFC", "p_val")
degenes$gsym <- row.names(degenes)
head(degenes)

### 2 plot
# 设置参数
logFCcut <- 0.56 #log2-foldchange
pvalCut <- 0.05 #P.value

#for dot size
logFCcut2 <- 1
logFCcut3 <- 2
pvalCut2 <- 0.01
pvalCut3 <- 0.001

### plot
n1 <- length(degenes[, 1])
cols <- rep("grey", n1)
names(cols)<- rownames(degenes)

#不同阈值的点的颜色
cols[degenes$p_val < pvalCut & degenes$avg_logFC > logFCcut]<- "#ca0000"
cols[degenes$p_val < pvalCut & degenes$avg_logFC < -logFCcut]<- "#2166ac"
table(cols)
color_transparent <- adjustcolor(cols, alpha.f = 0.5)
degenes$color_transparent <- color_transparent

# 複雜的的setting for size
n1 <- length(degenes[, 1])
size <- rep(1, n1)

#不同阈值的点的大小
size[degenes$p_val < pvalCut & degenes$avg_logFC > logFCcut]<- 2
size[degenes$p_val < pvalCut2 & degenes$avg_logFC > logFCcut2]<- 4
size[degenes$p_val < pvalCut3 & degenes$avg_logFC > logFCcut3]<- 6
size[degenes$p_val < pvalCut & degenes$avg_logFC < -logFCcut]<- 2
size[degenes$p_val < pvalCut2 & degenes$avg_logFC < -logFCcut2]<- 4
size[degenes$p_val < pvalCut3 & degenes$avg_logFC < -logFCcut3]<- 6
table(size)

# xlim
xmin <- min(degenes$avg_logFC)
xmax <- max(degenes$avg_logFC)
degenes$logp <- -log10(degenes$p_val)
degenes$logp[degenes$logp > 303] <- 303
# Construct the plot object
p1 <- ggplot(data=degenes, aes(avg_logFC, logp, label = gsym)) +
  geom_point(alpha = 0.6, size = size, colour = degenes$color_transparent) +
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  # xlim(c(ymin,ymax)) + 
  #画阈值分界线
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + #虚线的形状和粗细
  geom_hline(yintercept = -log10(pvalCut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_rect_size = 1.1,  # 修改坐标轴的外观，变成四方形的框，
           base_line_size = 1.1) + # 修改主题中默认的边框的粗细
  theme(panel.grid =element_blank(), #去除背景
        panel.border = element_rect(color = "black"), #之前边框的color是灰色的,这里调成黑色
        axis.text = element_text(size = 16, color = "black"), #调整坐标轴字体大小
        axis.title = element_text(size = 18, color = "black")) 

p1

### 添加文字
n = 1
p1 + geom_text_repel(aes(x = avg_logFC, y = -log10(p_val), 
                         label = ifelse(avg_logFC > n & p_val < 0.05, rownames(degenes),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines")) +
  geom_text_repel(aes(x = avg_logFC, y = -log10(p_val), 
                      label = ifelse(avg_logFC < -n & p_val < 0.05, rownames(degenes),"")),
                  colour="darkblue", size = 5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines"))


ggsave(filename = "Fig2__P12345_Primary_vs_metastasis_volcano.pdf", width = 6, height = 5.8)
#-----





#' 2021/07/18
#' DEgenes
#' P12345
#' CTC/P + CTC/M + P/M
#' top50
#' heatmap
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

dirs <- "~/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)

seuObject <- readRDS(file = "../P12345_pc50_res1_plot.rds")
table(seuObject$label)
Idents(seuObject) <- "label"

### 2 DE genes
degene1 <- read.csv(file = "../../DEgenes/P_M/P12345_eip_Primary_vs_Metastasis_degenes_0.5.csv", 
                    header = T, row.names = 1)
degene1 <- degene1[order(degene1$avg_log2FC, decreasing = T), ]

degene2 <- read.csv(file = "../../DEgenes/CTC/20210718_CTC_vs_Primary_DEgenes.csv", 
                    header = T, row.names = 1)
degene3 <- read.csv(file = "../../DEgenes/CTC/20210718_CTC_vs_Metastasis_DEgenes.csv", 
                    header = T, row.names = 1)

library(ggsci)
library(ComplexHeatmap)
library(circlize)

cells <- c("CTC", "Epithelial")
seu.sub <- subset(seuObject, subset = label %in% cells)

### 1 Plot : heatmap
group <- as.character(seu.sub$group_tissue)
table(group)
group[group=="Portal-Venous"] <- "CTC"

# plot : heatmap
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

order1 <- c("CTC", "Primary", "Metastasis")
orders <- c("CTC_P", "CTC_M", "Primary")
# colors
cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

mypal.tissue <- c("#D14A3D", "#62BFBF", "#62A1D9")
names(mypal.tissue) <- order1

# annotation
ha.t <- HeatmapAnnotation(Group = factor(group, levels = order1),
                          col = list(Group = mypal.tissue))

# genes
top = 50
degene1$gene <- rownames(degene1)
degene2$gene <- rownames(degene2)
degene3$gene <- rownames(degene3)

topn <- rbind(degene2[1:top, ], degene3[1:top, ], degene1[1:top, ])
topn$cluster <- rep(orders, c(50, 50, 50))
topn <- topn[!duplicated(topn$gene), ]

# data to plot
ht_opt$message = FALSE
count <- as.matrix(GetAssayData(object = seu.sub, slot = "data"))[rownames(topn), , drop = FALSE]
exMatrix <- t(scale(t(count)))

# plot
width = 10
height <- nrow(topn) * 0.1 + 2
pdf("Fig2_P12345_CTC_P_M_Degenes.pdf", height = height, width = width)
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
        column_title_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        row_names_max_width = max_text_width(rownames(exMatrix),
                                             gp = gpar(fontsize = 10)))
dev.off()
#-----




#' 2021/09/07
#' DEgenes
#' P12345
#' CTC/P+M epi
#' top30
#' heatmap
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

dirs <- "~/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)

seuObject <- readRDS(file = "../P12345_pc50_res1_plot.rds")
table(seuObject$label)
Idents(seuObject) <- "label"

### 2 DE genes
degene1 <- read.csv(file = "../../DEgenes/CTC/20210811_CTC_vs_epithelial_DEgenes.csv", 
                    header = T, row.names = 1)
degene1 <- degene1[order(degene1$avg_log2FC, decreasing = T), ]

library(ggsci)
library(ComplexHeatmap)
library(circlize)

cells <- c("CTC", "Epithelial")
seu.sub <- subset(seuObject, subset = label %in% cells)

### 1 Plot : heatmap
group <- as.character(seu.sub$group_tissue)
table(group)
group[group=="Portal-Venous"] <- "CTC"

# plot : heatmap
library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

order1 <- c("CTC", "Primary", "Metastasis")
orders <- c("CTC", "Epithelial")
# colors
cols = colorRamp2(seq(from = -2, to = 2, length.out = 11),
                  c("#053061", "#2166AC", "#4393C3",
                    "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7",
                    "#F4A582", "#D6604D","#B2182B", "#67001F"))

mypal.tissue <- c("#D14A3D", "#62BFBF", "#62A1D9")
names(mypal.tissue) <- order1

# annotation
ha.t <- HeatmapAnnotation(Group = factor(group, levels = order1),
                          col = list(Group = mypal.tissue))

# genes
top = 30
degene1$gene <- rownames(degene1)
topn <- degene1[c(1:top, (nrow(degene1)-(top-1)):nrow(degene1)), ]
topn$cluster <- rep(orders, c(top, top))

# data to plot
ht_opt$message = FALSE
count <- as.matrix(GetAssayData(object = seu.sub, slot = "data"))[rownames(topn), , drop = FALSE]
exMatrix <- t(scale(t(count)))

# plot
width = 9
height <- nrow(topn) * 0.2 + 2
pdf("Fig2_P12345_CTC_vs_epithelial_Degenes_top30.pdf", height = height, width = width)
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




#' 2022/01/20
#' enrichment analysis
#' CTC
#' selected genesets of CTC
#' Heatmap
#-----
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

dirs <- "~/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)

seuObject <- readRDS(file = "../P12345_pc50_res1_plot.rds")

library(ggplot2)
library(ggsci)
library(ComplexHeatmap)
library(circlize)

tissue <- "Pancreatic_P12345"
data <- "20210719"

## 1 subset
cells <- c("Epithelial", "CTC")
seu.sub <- subset(seuObject, subset = label %in% cells)

## 2 input gsva score
# 使用fread读入的表格是没有行名的,需要先把第一列添加为行名,然后再删除第一列.
score1 <- fread(file = paste0("../../scss_output/H_", tissue,"_data.tsv"), header = T, data.table = F)
score2 <- fread(file = paste0("../../scss_output/C2_CP_KEGG_", tissue,"_data.tsv"), header = T, data.table = F)
score3 <- fread(file = paste0("../../scss_output/C5_BP_", tissue,"_data.tsv"), header = T, data.table = F)
score4 <- fread(file = paste0("../../scss_output/C6_", tissue,"_data.tsv"), header = T, data.table = F)
score5 <- fread(file = paste0("../../scss_output/C2_CP_REACTOME_", tissue,"_data.tsv"), header = T, data.table = F)
gsvascore <- cbind(score1, score2, score3, score4, score5)
gsvascore[1:5, 1:5]

## 3 geneset to plot
genesets <- read.csv(file = "../../data/CTC_genesets_heatmap_20220210.csv", header = T)
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

write.csv(df.plot, file = "20220210_Fig1_heatmap_genesets.csv")
write.csv(group, file = "20220210_Fig1_heatmap_genesets_group.csv")

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
pdf("20220120_Fig1_P12345_CTC_pathays_heatmap.pdf", height = height, width = width)
Heatmap(df.plot, 
        cluster_rows = F, 
        cluster_columns = F)
lgd = Legend(labels = group.row$Var1, title = "Genesets", legend_gp = gpar(fill = col))
draw(lgd, x = unit(0.9, "npc"), y = unit(0.9, "npc"), just = c("right", "top"))
dev.off()


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
#-----



#' 2021/09/20
#' heatmap
#' CTC
#' 说明CTC具有血小板基因
#-----
dirs <- "/home/SJE/project/pancreatic_scRNA/analysis/Figures/Figure2/CTC_gene"
setwd(dirs)

genes <- c('PPBP','B2M','FTL','PF4','SH3BGRL3','RPS27','PTMA','OAZ1', "ITGA2B", "ITGB3", "RGS18")

seuObject <- readRDS("~/project/pancreatic_scRNA/analysis/Figures/P12345_non_immune_pc40_res0.5.rds")
seu.ctc <- subset(seuObject, subset = label2 %in% c("CTC"))
seu.ctc <- subset(seu.ctc, subset = PTPRC = 0)

df <- data.frame(gene=as.matrix(GetAssayData(seu.ctc, slot = "data"))[genes, ])
library(ComplexHeatmap)

# colors
cols = colorRamp2(seq(from = 0, to = 5, length.out = 4),
                  c("white", "#FDDBC7","#F4A582", "red"))

pdf(file = "Fig2S_CTC_platelet_genes_heatmap.pdf", width = 15, height = 4)
Heatmap(df, show_column_names = F, col = cols, rect_gp = gpar(col= "lightgrey"),
        name = "Log2(TPM)")
dev.off()
#-----


#' 2021/09/08
#' boxplot
#' CTC+ primary + metastasis
#' 说明CTC具有血小板基因和上皮基因的丢失以及其它特征
#-----
dirs <- "/home/SJE/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)
dir.create("CTC_gene")
library(ggplot2)
library(ggpubr)

genes <- c("PPBP", "PF4", "GP9", "ITGA2B",
           'NRGN','RGS18', 'CCL5', 'GNG11', 
           'SRGN', 'TGFB1', 'SH3BGRL3', 'SPARC', 
           "CDH1", "KRT8", "KRT18", "KRT19")
p12345 <- readRDS(file = "../P12345_pc50_res1_plot.rds")
seuObject <- subset(p12345, subset = label2 %in% c("CTC", "Epithelial"))
group <- as.character(seuObject$group_tissue)
table(group)
group[group == "Portal-Venous"] <- "CTC"

orders <- c("Primary","CTC", "Metastasis")
Idents(seuObject) <- factor(group, levels = orders)

# df to plot
df.plot <- as.data.frame(t(as.matrix(GetAssayData(seuObject))[genes, ]))
head(df.plot)

df.plot$group <- Idents(seuObject)
# out
write.csv(df.plot, file = "20220114_FigS5d_boxplot.csv")

## plot
lplot <- list()
for (i in 1:length(genes)) {
  name <- genes[i]
  p1 <- ggboxplot(df.plot, x = "group", y = name, fill = "group", outlier.shape = NA,
                  palette = c("#62BFBF", "#ca0000", "#62A1D9"), bxp.errorbar = T) +
    xlab("") + ylab("log2(TPM + 1)") +
    ggtitle(label = name) +
    guides(fill = F) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  lplot[[i]] <- p1
  
}

lplot[[11]]

# plot
lplot[[1]]+lplot[[2]]+lplot[[3]]+lplot[[4]]+lplot[[5]]+lplot[[6]]+lplot[[7]]+lplot[[8]]+
  lplot[[9]]+lplot[[10]]+lplot[[11]]+lplot[[12]]+lplot[[13]]+lplot[[14]]+lplot[[15]]+
  lplot[[16]]+plot_layout(nrow = 4)
ggsave(filename = paste0("CTC_gene/Fig2_P12345_CTC_genes_boxplot.pdf"), width = 5.2, height = 12)
#-----



#' 2021/09/08
#' 说明CTC中有血小板是普遍的现象
#' 找到数据集：
#' (1)https://ftp.cngb.org/pub/CNSA/data3/CNP0000095/Single_Cell/CSE0000014/
#' (2)https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144561
#' (3)https://www.sciencedirect.com/science/article/pii/S2211124714007050  小鼠数据
#-----
dirs <- "/home/SJE/project/pancreatic_scRNA/analysis/Figures/Figure2/"
setwd(dirs)

library(Seurat)
library(ggpubr)

genes <- c('HLA-E','PPBP','B2M','FTL','PF4','SH3BGRL3','RPS27','PTMA','OAZ1')

### 1 肝癌CTC单细胞数据 NC 2021
df.ctc <- read.csv(file = "../../data/CTC/CTC_HCC_log_tpm_expression_matrix.txt", sep = "\t", 
                   header = T, row.names = 1)
df.ctc <- df.ctc[genes, ]

group <- unlist(lapply(strsplit(colnames(df.ctc), split="[.]"), function(x) x[2]))

# plot
df.plot <- as.data.frame(t(df.ctc))
tail(df.plot)

df.plot$group <- group

lplot <- list()
for (i in 1:length(genes)) {
  name <- genes[i]
  p1<-ggboxplot(df.plot, x = "group", y = name, fill = "group", bxp.errorbar = T, 
                palette = pal_npg()(4)) +
    xlab("") + ylab("log2(TPM)") +
    ggtitle(label = name) +
    guides(fill = F) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))
  lplot[[i]] <- p1
}
lplot[[1]]+lplot[[2]]+lplot[[3]]+lplot[[4]]+lplot[[5]]+lplot[[6]]+lplot[[7]]+lplot[[8]]+
  plot_layout(nrow = 1)
ggsave(filename = paste0("CTC_gene/Fig2_HCC_CTC_Platelet_genes.pdf"), width = 12, height = 3)


### 2 PAAD CTC
df.ctc <- read.csv(file = "../../data/CTC/PAAD_CTC_GSE144561_rawCountsAllsamples.txt", sep = "\t", 
                   header = T, row.names = 1)
colnames(df.ctc)
# add group info
df.ctc <- df.ctc[, 1:56]
group <- rep(c("HD", "locPDAC", "metPDAC"), c(21, 17, 18))
df.ctc <- t(t(df.ctc)*10000/colSums(df.ctc))
df.ctc <- df.ctc[genes, ]
df.ctc <- log2(df.ctc+1)

# plot
df.plot <- as.data.frame(t(df.ctc))
tail(df.plot)

df.plot$group <- group

# plot
lplot <- list()
for (i in 1:length(genes)) {
  name <- genes[i]
  p1<-ggboxplot(df.plot, x = "group", y = name, fill = "group", bxp.errorbar = T, 
                palette = pal_npg()(4)) +
    xlab("") + ylab("log2(TPM)") +
    ggtitle(label = name) +
    guides(fill = F) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))
  lplot[[i]] <- p1
}
lplot[[1]]+lplot[[2]]+lplot[[3]]+lplot[[4]]+lplot[[5]]+lplot[[6]]+lplot[[7]]+lplot[[8]]+
  plot_layout(nrow = 1)
ggsave(filename = paste0("CTC_gene/Fig2_PDAC_CTC_Platelet_genes.pdf"), width = 12, height = 3)




### 3 mouse CTC
df.ctc <- read.csv(file = "../../data/CTC/PAAD_mouse_GSE51372_count_genes_filter_20210104.csv", 
                   header = T, row.names = 1)
# groups
df.group <-read.csv(file = "../../data/CTC/PAAD_mouse_GSE51372_groups.csv", header = T)
df.group$sample <- gsub("-", ".", df.group$sample)
table(df.group$groups)
df.ctc <- df.ctc[, df.group$groups %in% c("CTC_c", "CTC_plt", "Tumor")]
df.group <- df.group[df.group$groups %in% c("CTC_c", "CTC_plt", "Tumor"), ]

df.ctc <- t(t(df.ctc)*10000/colSums(df.ctc))
genes <- c("Ppbp", "B2m", "Ftl1","Pf4","Sh3bgrl3","Rps27","Ptma" ,"Oaz1")
df.ctc <- df.ctc[genes, ]
df.ctc <- log2(df.ctc+1)

# plot
df.plot <- as.data.frame(t(df.ctc))
tail(df.plot)
group <- df.group$groups

df.plot$group <- factor(group, levels = c("CTC_c", "CTC_plt", "Tumor"))

# plot
lplot <- list()
for (i in 1:length(genes)) {
  name <- genes[i]
  p1<-ggboxplot(df.plot, x = "group", y = name, fill = "group", bxp.errorbar = T, 
                palette = pal_npg()(4)) +
    xlab("") + ylab("log2(TPM)") +
    ggtitle(label = name) +
    guides(fill = F) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))
  lplot[[i]] <- p1
}
lplot[[1]]+lplot[[2]]+lplot[[3]]+lplot[[4]]+lplot[[5]]+lplot[[6]]+lplot[[7]]+lplot[[8]]+
  plot_layout(nrow = 1)
ggsave(filename = paste0("CTC_gene/Fig2_mouse_PDAC_CTC_Platelet_genes.pdf"), width = 12, height = 3)
#-----


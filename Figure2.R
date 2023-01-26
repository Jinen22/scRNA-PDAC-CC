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



#' Figure 2
#' Figure 2D
#' cellphoneDB
#' network
#' igraph
#' Pirmary
#-----
library(psych)
library(igraph)
library(tidyverse)
library(ggsci)

## colors
cellorder <- c("B cell", "Naive T", "CD8 EFF", "CD8 Ex", "Memory T", "Treg", "NKT", "NK", 
               "cDC", "pDC", "M1", "M2", "Monocyte", "Neutrophil", "Mast", "CTC", "Primary", 
               "Metastasis")
# celltype
length(cellorder)
mypal_cell <- c(pal_npg()(7), 
                c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3",
                  "#F781BE", "#33A9CF", "#E51A1D", "#62BFBF", "#62A1D9"))
names(mypal_cell) <- cellorder

tissue <- "Primary"


### 1 Primary
files = "./Primary/"
df.p = read.table(paste0(files, "/outs/pvalues.txt"), header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
df.m <- read.table(paste0(files, "/outs/means.txt"), header=T, stringsAsFactors = F, sep = '\t', 
                   comment.char = '', check.names=F)
rownames(df.p) <- df.p$interacting_pair
rownames(df.m) <- df.m$interacting_pair

### 2 data filtered
df.p <- df.p[, -(1:11)]
df.m <- df.m[, -(1:11)]

# keep connect rows
df <- read.csv(file = "../../Figures/Figure3/igraph/Cellchat_cell_contact_gene_20211213.csv", header = T)


keep <- rownames(df.p) %in% df$interaction_name
table(keep)
df.p <- df.p[keep, ]
df.m <- df.m[keep, ]

df.plot <- c()
for (cell in colnames(df.p)) {
  pval <- df.p[, cell]
  mea <- df.m[, cell]
  mea <- mea[pval < 0.05]
  df.plot <- c(df.plot, sum(mea))
}
names(df.plot) <- colnames(df.p)
df.plot

name1 <- unlist(lapply(strsplit(names(df.plot), split = "[|]"), function (x) x[1]))
name2 <- unlist(lapply(strsplit(names(df.plot), split = "[|]"), function (x) x[2]))
df.plot <- data.frame(source  = name1, target = name2, count = df.plot)
head(df.plot)
keep <- df.plot$source != df.plot$target
df.plot <- df.plot[keep, ]
df.plot$count <- round(df.plot$count, digits = 2)

# plot
mynet <- df.plot
table(mynet$count)
mynet %>% filter(count>0) -> mynet  # filter
head(mynet)
net<- graph_from_data_frame(mynet) 

### 3 plot
length(unique(mynet$source))
## CTC
# colors
df.order <- mynet[grep("^CTC", mynet$source), ]
df.order <- df.order[order(df.order$count, decreasing = T), ]
grouporder <- c("CTC", df.order$target)
allcolour=mypal_cell[grouporder]
allcolour

coords <- layout_in_circle(net, order = grouporder)  
E(net)$width  <- E(net)$count/1  

#plot
length(unique(mynet$source))
net1<-net

i <- tissue
E(net1)$count <- ""
E(net1)[map(grouporder,function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$count  <- E(net)[map(grouporder,function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$count  

E(net1)[map(unique(mynet$source),function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$color <- allcolour[i]

pdf(file = "Primary_CTC_cellphoneDB_cell_contect_network.pdf", width = 5, height = 5)
plot(net1, edge.arrow.size=.1, 
     edge.curved=0.2, 
     edge.label = E(net1)$count, 
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=1.5,
     vertex.label.dist=3
) 
dev.off()
#-----



#' Figure 2
#' Figure 2E
#' cellphoneDB
#' network
#' igraph
#' Blood
#-----
library(psych)
library(igraph)
library(tidyverse)
library(ggsci)

## colors
cellorder <- c("B cell", "Naive T", "CD8 EFF", "CD8 Ex", "Memory T", "Treg", "NKT", "NK", 
               "cDC", "pDC", "M1", "M2", "Monocyte", "Neutrophil", "Mast", "CTC", "Primary", 
               "Metastasis")
# celltype
length(cellorder)
mypal_cell <- c(pal_npg()(7), 
                c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3",
                  "#F781BE", "#33A9CF", "#E51A1D", "#62BFBF", "#62A1D9"))
names(mypal_cell) <- cellorder


### 1 blood
files = "./Blood"
df.p = read.table(paste0(files, "/outs/pvalues.txt"), header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
df.m <- read.table(paste0(files, "/outs/means.txt"), header=T, stringsAsFactors = F, sep = '\t', 
                   comment.char = '', check.names=F)
rownames(df.p) <- df.p$interacting_pair
rownames(df.m) <- df.m$interacting_pair

### 2 data filtered
df.p <- df.p[, -(1:11)]
df.m <- df.m[, -(1:11)]

# keep connect pairs
df <- read.csv(file = "../../Figures/Figure3/igraph/Cellchat_cell_contact_gene_20211213.csv", header = T)
keep <- rownames(df.p) %in% df$interaction_name
table(keep)
df.p <- df.p[keep, ]
df.m <- df.m[keep, ]


df.plot <- c()
for (cell in colnames(df.p)) {
  pval <- df.p[, cell]
  mea <- df.m[, cell]
  mea <- mea[pval < 0.05]
  df.plot <- c(df.plot, sum(mea))
}
names(df.plot) <- colnames(df.p)
df.plot

name1 <- unlist(lapply(strsplit(names(df.plot), split = "[|]"), function (x) x[1]))
name2 <- unlist(lapply(strsplit(names(df.plot), split = "[|]"), function (x) x[2]))
df.plot <- data.frame(source  = name1, target = name2, count = df.plot)
head(df.plot)
keep <- df.plot$source != df.plot$target
df.plot <- df.plot[keep, ]
df.plot$count <- round(df.plot$count, digits = 2)

# plot
mynet <- df.plot
table(mynet$count)
mynet %>% filter(count>0) -> mynet  # filter
head(mynet)
net<- graph_from_data_frame(mynet) 

### 3 plot
length(unique(mynet$source))
## CTC
# colors
df.order <- mynet[grep("^CTC", mynet$source), ]
df.order <- df.order[order(df.order$count, decreasing = T), ]
grouporder <- c("CTC", df.order$target)
allcolour=mypal_cell[grouporder]
allcolour

coords <- layout_in_circle(net, order = grouporder)  
E(net)$width  <- E(net)$count/1  

#plot
net1<-net

i <- "CTC"
E(net1)$count <- ""
E(net1)[map(grouporder,function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$count  <- E(net)[map(grouporder,function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$count  

E(net1)[map(unique(mynet$source),function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$color <- allcolour[i]

pdf(file = "Blood_CTC_cellphoneDB_cell_contect_network.pdf", width = 5, height = 5)
plot(net1, edge.arrow.size=.1, 
     edge.curved=0.2, 
     edge.label = E(net1)$count,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=1.5,
     vertex.label.dist=3
) 
dev.off()
#-----







#' Figure 2
#' Figure 2F
#' cellphoneDB
#' network
#' igraph
#' Metastasis
#-----
library(psych)
library(igraph)
library(tidyverse)
library(ggsci)

## colors
cellorder <- c("B cell", "Naive T", "CD8 EFF", "CD8 Ex", "Memory T", "Treg", "NKT", "NK", 
               "cDC", "pDC", "M1", "M2", "Monocyte", "Neutrophil", "Mast", "CTC", "Primary", 
               "Metastasis")
# celltype
length(cellorder)
mypal_cell <- c(pal_npg()(7), 
                c("#377DB8", "#4DAF4A", "#A4572B", "#E51A1D", "#FF7F00", "#984EA3",
                  "#F781BE", "#33A9CF", "#E51A1D", "#62BFBF", "#62A1D9"))
names(mypal_cell) <- cellorder

tissue <- "Metastasis"


### 1 Metastasis
files = "./Metastasis/"
df.p = read.table(paste0(files, "/outs/pvalues.txt"), header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
df.m <- read.table(paste0(files, "/outs/means.txt"), header=T, stringsAsFactors = F, sep = '\t', 
                   comment.char = '', check.names=F)
rownames(df.p) <- df.p$interacting_pair
rownames(df.m) <- df.m$interacting_pair

### 2 data filtered
df.p <- df.p[, -(1:11)]
df.m <- df.m[, -(1:11)]

# keep connect pairs
df <- read.csv(file = "../../Figures/Figure3/igraph/Cellchat_cell_contact_gene_20211213.csv", header = T)


keep <- rownames(df.p) %in% df$interaction_name
table(keep)
df.p <- df.p[keep, ]
df.m <- df.m[keep, ]

df.plot <- c()
for (cell in colnames(df.p)) {
  pval <- df.p[, cell]
  mea <- df.m[, cell]
  mea <- mea[pval < 0.05]
  df.plot <- c(df.plot, sum(mea))
}
names(df.plot) <- colnames(df.p)
df.plot

name1 <- unlist(lapply(strsplit(names(df.plot), split = "[|]"), function (x) x[1]))
name2 <- unlist(lapply(strsplit(names(df.plot), split = "[|]"), function (x) x[2]))
df.plot <- data.frame(source  = name1, target = name2, count = df.plot)
head(df.plot)
keep <- df.plot$source != df.plot$target
df.plot <- df.plot[keep, ]
df.plot$count <- round(df.plot$count, digits = 2)


### plot
mynet <- df.plot

table(mynet$count)
mynet %>% filter(count>0) -> mynet 
head(mynet)
net<- graph_from_data_frame(mynet) 

### 3 Network
# colors
df.order <- mynet[grep(tissue, mynet$source), ]
df.order <- df.order[order(df.order$count, decreasing = T), ]
grouporder <- c(tissue, df.order$target)
allcolour=mypal_cell[grouporder]

coords <- layout_in_circle(net, order = grouporder)  
E(net)$width  <- E(net)$count/1  

#plot
length(unique(mynet$source))
net1<-net

i <- tissue
E(net1)$count <- ""
E(net1)[map(grouporder,function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$count  <- E(net)[map(grouporder,function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$count  

E(net1)[map(unique(mynet$source),function(x) {
  get.edge.ids(net,vp = c(i,x))
})%>% unlist()]$color <- allcolour[i]

pdf(file = "Metastasis_CTC_cellphoneDB_cell_contect_network.pdf", width = 5, height = 5)
plot(net1, edge.arrow.size=.1, 
     edge.curved=0.2, 
     edge.label = E(net1)$count, 
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=1.5,
     vertex.label.dist=3
) 
dev.off()
#-----





#' Figure 2
#' Figure 2G
#' cellphoneDB
#' (1) immune check point genes
#' (2) chemokine genes
#' Primary + CTC + Metastasis
#' plot
#-----
library(Seurat)
library(data.table)

df <- read.csv(file = "../data/cellphonedb_gene.csv", header = T)
unique(df$genes)

### 1 import cellphoneDB results
cellorder <- c("B cell", "Naive T", "CD8 EFF", "CD8 Ex", "Memory T", "Treg", "NKT", "NK", 
               "cDC", "pDC", "M1", "M2", "Monocyte", "Neutrophil", "Mast")
cellread <- function(files, tissue, cellorder) {
  df <- read.table(files, header=T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names=F)
  df <- df[, -c(1, 3:11)]
  df <- df[, c(1, grep(tissue, colnames(df)))]
  colnames(df)
  filter <- which(colnames(df) == paste0(tissue, "|", tissue))
  df <- df[, -filter]
  # order
  colorder <- c(colnames(df)[1], 
                paste0(tissue, "|", cellorder), paste0(cellorder, "|", tissue))
  colorder <- colorder[colorder %in% colnames(df)]
  df <- df[, colorder]
  
  return(df)
}

# pvalue
pvalue.p <- cellread("cellphonedb/Primary/outs/pvalues.txt", tissue = "Primary", 
                     cellorder = cellorder)
pvalue.m <- cellread("cellphonedb/Metastasis/outs/pvalues.txt", tissue="Metastasis", 
                     cellorder = cellorder)
pvalue.b <- cellread("cellphonedb/Blood/outs/pvalues.txt", tissue = "CTC", 
                     cellorder = cellorder)
df.p <- full_join(pvalue.b, pvalue.p, by="interacting_pair")
df.p <- full_join(df.p, pvalue.m, by="interacting_pair")
df.p[is.na(df.p)] <- 1

# mean
mean.p <- cellread("cellphonedb/Primary/outs/means.txt", tissue = "Primary", 
                   cellorder = cellorder)
mean.m <- cellread("cellphonedb/Metastasis/outs/means.txt", tissue="Metastasis", 
                   cellorder = cellorder)
mean.b <- cellread("cellphonedb/Blood/outs/means.txt", tissue = "CTC", 
                   cellorder = cellorder)
df.m <- full_join(mean.b, mean.p, by="interacting_pair")
df.m <- full_join(df.m, mean.m, by="interacting_pair")
df.m[is.na(df.m)] <- 0

#### 1 immune genes
gene.group <- "immune"
genes <- unique(df$genes[df$group == gene.group])
#1 rows filter
keep <- c()
for (i in genes) {
  pair <- grep(i, df.p$interacting_pair)
  keep <- c(keep, pair)
}
keep.row <- unique(keep)
pairs <- df.p$interacting_pair[keep.row]

# 2 pvalue
df.p.use <- df.p[keep.row, -1]
rownames(df.p.use) <- pairs
df.m.use <- df.m[keep.row, -1]
rownames(df.m.use) <- pairs

### 2 plot
df.pvalue.sub <- df.p.use
df.mean.sub <- df.m.use
# 2.1 filter by pvalue
row.filter1 <- -which(apply(df.pvalue.sub, 1, function(x) all(x > 0.05)))
df.pvalue.use <- df.pvalue.sub[row.filter1, ]
df.mean.use <- df.mean.sub[row.filter1, ]

## 2.2 filter by mean
row.filter2 <- -which(apply(df.mean.use, 1, function(x) all(x < 0.5)))
df.pvalue.use <- df.pvalue.use[row.filter2, ]
df.mean.use <- df.mean.use[row.filter2, ]

## 3.1 rownames and colnames
pair <- rownames(df.pvalue.use)
cluster <- colnames(df.pvalue.use)

df.use <- expand.grid(pair, cluster)
head(df.use)

## 3.2 pvalue and mean
pval <- unlist(df.pvalue.use)
pval[pval==0] <- 0.0009
plot.data <- cbind(df.use, -log10(pval))

pr <- unlist(df.mean.use)
# pr[pr > 2] <- 2

plot.data = cbind(plot.data, pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
head(plot.data)

# add group
plot.data$group <- rep(c("CTC", "Primary", "Metastasis"), 
                       c(12*length(pair), 28*length(pair), 28*length(pair)))
plot.data$group <- factor(plot.data$group, levels = c("Primary", "CTC", "Metastasis"))
## plot

# 配色
my_palette <- my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
width=0.23
witdh <- width * length(colnames(df.pvalue.sub)) + 4
height <-  0.2*length(pair) + 3

pdf(file = "Immune_checkpoint_cellphoneDB.pdf", width = witdh, height = height)
ggplot(plot.data, aes(x=clusters, y=pair)) +
  geom_point(aes(size=pvalue, color=mean)) +
  scale_color_gradientn('Mean (Molecule 1, Molecule 2)', 
                        colors=my_palette) + # colors
  theme_bw() +
  facet_grid(.~group, scales = "free", space = "free") +
  theme(panel.grid.minor = element_blank(), # blank
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

dev.off()


#### 2 chemokine genes
gene.group <- "chemokine"
genes <- unique(df$genes[df$group == gene.group])
#1 rows filter
keep <- c()
for (i in genes) {
  pair <- grep(i, df.p$interacting_pair)
  keep <- c(keep, pair)
}
keep.row <- unique(keep)
pairs <- df.p$interacting_pair[keep.row]

# 2 pvalue
df.p.use <- df.p[keep.row, -1]
rownames(df.p.use) <- pairs
df.m.use <- df.m[keep.row, -1]
rownames(df.m.use) <- pairs

### 2 plot
# pairs : interection pair names of df.pvalue.sub and df.mena.sub
df.pvalue.use <- df.p.use
df.mean.use <- df.m.use

## 3.1 rownames and colnames
pair <- rownames(df.pvalue.use)
cluster <- colnames(df.pvalue.use)

df.use <- expand.grid(pair, cluster)
head(df.use)

## 3.1 pvalue and mean
pval <- unlist(df.pvalue.use)
pval[pval==0] <- 0.0009
plot.data <- cbind(df.use, -log10(pval))

pr <- unlist(df.mean.use)
# pr[pr > 2] <- 2

plot.data = cbind(plot.data, pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
head(plot.data)

# add group
plot.data$group <- rep(c("CTC", "Primary", "Metastasis"), 
                       c(12*length(pair), 28*length(pair), 28*length(pair)))
plot.data$group <- factor(plot.data$group, levels = c("Primary", "CTC", "Metastasis"))
## plot

# color
my_palette <- my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
width=0.23
witdh <- width * length(colnames(df.pvalue.sub)) + 4
height <-  0.2*length(pair) + 3

pdf(file = "20220713_chemokine_cellphoneDB.pdf", width = witdh, height = height)
ggplot(plot.data, aes(x=clusters, y=pair)) +
  geom_point(aes(size=pvalue, color=mean)) +
  scale_color_gradientn('Mean (Molecule 1, Molecule 2)', 
                        colors=my_palette) + # 用自定义颜色画点
  theme_bw() +
  facet_grid(.~group, scales = "free", space = "free") +
  theme(panel.grid.minor = element_blank(), #不画网格
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) #边框

dev.off()
#-----




#' Figure 2
#' Figure 2H and 2J
#' barplot
#-----
library(Seurat)
library(dplyr)
library(ggpubr)
library(patchwork)

genes <- c('HLA-E', 'KLRD1')

seuOject <- readRDS(file = "Pall_tumor_ctc_immune.rds")
table(seuOject$label)

### 1 HLA-E in scRNA data
# epi + ctc
seu.sub <- subset(seuOject, subset = label %in% c("Tumor cells", "CTC"))
group1 <- seu.sub$group_tissue
table(group1)
group1[group1 == "Portal-Venous"] <- "CTC"
group1[group1 == "Primary"] <- "Primary tumor"
group1[group1 == "Metastasis"] <- "Metastasis tumor"
Idents(seu.sub) <- factor(group1, levels = c("Primary tumor", "CTC", "Metastasis tumor"))

df.plot <- as.data.frame(t(as.matrix(GetAssayData(seu.sub, slot = "data")[genes, ])))
head(df.plot)
df.plot$group <- factor(group1, levels = c("Primary tumor", "CTC", "Metastasis tumor"))

## plot
ggbarplot(df.plot, x = "group", y=genes[1], add="mean_se",fill="group",  position=position_dodge(0.8),
          palette = c("#62BFBF", "#ca0000", "#62A1D9")) +
  xlab("Cell type") +
  ylab("Expression") +
  ggtitle(genes[1]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(fill = F)

ggsave("Fig2_HLA-E_CTC_P_M_barplot.pdf", width = 2.5, height = 4)


### 2 GEO data
df <- as.matrix(GetAssayData(seu.sub, slot = "data")[genes, ])

# 1 HCC CTC scRNA
df.ctc1 <- read.csv(file = "../data/CTC/CTC_HCC_log_tpm_expression_matrix.txt", sep = "\t", 
                    header = T, row.names = 1)
df.ctc1 <- df.ctc1[genes, ]
df.ctc1 <- log2(df.ctc1 + 1)

# PDAC CTC
df.ctc2 <- read.csv(file = "../data/CTC/PAAD_CTC_GSE144561_rawCountsAllsamples.txt", sep = "\t", 
                    header = T, row.names = 1)
# remove healthy donor
filter <- grep("^HD", colnames(df.ctc2))
df.ctc2 <- df.ctc2[, -filter]
df.ctc2 <- t(t(df.ctc2)*10000/colSums(df.ctc2))
df.ctc2 <- df.ctc2[genes, ]
df.ctc2 <- log2(df.ctc2+1)

#3 ctcRbase数据
myread <- function(files) {
  df.test <- read.csv(file = files, header = T, row.names = 1)
  df.test <- t(t(df.test)*10000/colSums(df.test))
  gene <- unlist(lapply(strsplit(rownames(df.test), split = "_"), function(x) x[2]))
  keep <- c(grep("HLA-E", gene), grep("HLA-C", gene))
  df.test <- df.test[keep, ]
  rownames(df.test) <- genes
  df.test <- log2(df.test+1)
  return(df.test)
}
df.ctc3 <- myread(files = "../data/CTC/ctcRbase/LIHC_GSE117623_gene.mat.csv")
df.ctc4 <- myread(files = "../data/CTC/ctcRbase/breast_ctc_GSE67939.csv")
df.ctc5 <- myread(files = "../data/CTC/ctcRbase/breast_ctc_GSE86978.csv")
df.ctc6 <- myread(files = "../data/CTC/ctcRbase/color_GSE74369_gene_FPKM_mat.csv")
df.ctc7 <- myread(files = "../data/CTC/ctcRbase/mela_GSE38495_gene_FPKM_mat.csv")

genename <- unlist(lapply(strsplit(rownames(df.ctc7), split = "_"), function (x) x[2]))
df.test <- t(t(df.ctc7)*10000/colSums(df.ctc7))
gene <- unlist(lapply(strsplit(rownames(df.test), split = "_"), function(x) x[2]))
keep <- c(grep("HLA-E", gene), grep("HLA-C", gene))
df.test <- df.test[keep, ]
rownames(df.test) <- genes
df.test <- log2(df.test+1)
return(df.test)


# plot
df.plot <- as.data.frame(t(cbind(df,df.ctc1,df.ctc3,df.ctc2,df.ctc4,df.ctc5,df.ctc6)))
tail(df.plot)

group1[group1 == "CTC"] <- "PDAC CTC"
group1[group1 == "Primary tumor"] <- "PDAC primary tumor"
group1[group1 == "Metastasis tumor"] <- "PDAC metastasis tumor"
table(group1)

tissue <- c("CNP0000095 HCC CTC", "GSE117623 HCC CTC","GSE144561 PDAC CTC", "GSE67939 BRCA CTC",
            "GSE86978 BRCA CTC","GSE74369 COAD CTC")
df.plot$group <- c(group1, 
                   rep(tissue, 
                       c(ncol(df.ctc1), ncol(df.ctc3), ncol(df.ctc2),ncol(df.ctc4),ncol(df.ctc5),
                         ncol(df.ctc6))))
table(df.plot$group)
df.plot$group <- factor(df.plot$group, 
                        levels = c("PDAC primary tumor", "PDAC CTC", tissue,"PDAC metastasis tumor"))

# output
df.out <- df.plot[, c(1,3)]

## plot
my_cpmpare <- list(c("PDAC primary tumor", "PDAC CTC"), 
                   c("PDAC primary tumor", "CNP0000095 HCC CTC"),
                   c("PDAC primary tumor", "GSE117623 HCC CTC"),
                   c("PDAC primary tumor", "GSE144561 PDAC CTC"),
                   c("PDAC primary tumor", "GSE67939 BRCA CTC"),
                   c("PDAC primary tumor", "GSE86978 BRCA CTC"),
                   c("PDAC primary tumor", "GSE74369 COAD CTC"))

ggbarplot(df.plot, x = "group", y=genes[1], add="mean_se",fill="group",  position=position_dodge(0.8),
          palette = c("#62BFBF", rep("#ca0000", length(tissue)+1),"#62A1D9")) +
  xlab("Cell type") +
  ylab("Log2 (TPM + 1)") +
  ggtitle(genes[1]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(fill = F) +
  stat_compare_means(comparisons = my_cpmpare)
ggsave("Fig2_HLAE_GEO_CTC_barplot_pvalue.pdf", width = 5, height = 5)
#-----




#' Figure 2
#' Figure 2I
#' barplot
#' mouse CTC
#-----
df.ctc <- read.csv(file = "../../../data/CTC/PAAD_mouse_GSE51372_count_genes_filter_20210104.csv", 
                   header = T, row.names = 1)
# groups
df.group <-read.csv(file = "../../../data/CTC/PAAD_mouse_GSE51372_groups.csv", header = T)
df.group$sample <- gsub("-", ".", df.group$sample)
table(df.group$groups)
df.ctc <- df.ctc[, df.group$groups %in% c("CTC_c","CTC_plt","CTC_pro", "Tumor")]
df.group <- df.group[df.group$groups %in% c("CTC_c","CTC_plt","CTC_pro", "Tumor"), ]

df.ctc <- t(t(df.ctc)*1000000/colSums(df.ctc))
genes <- c("H2-T23", "Ppbp")
df.ctc <- df.ctc[genes, ]
df.ctc <- log2(df.ctc+1)

# plot
df.plot <- as.data.frame(t(df.ctc))
tail(df.plot)
group <- df.group$groups
group[group %in% c("CTC_c","CTC_plt","CTC_pro")] <- "CTC"
df.plot$group <- factor(group, levels = c("Tumor", "CTC"))

# output
df.out <- df.plot[, c(1,3)]

## plot
# P value
my_compare <- list(c("Tumor", "CTC"))
ggbarplot(df.plot, x = "group", y=genes[1], add="mean_se",fill="group",  position=position_dodge(0.8),
          palette = c("#62BFBF", "#ca0000", "#62A1D9")) +
  xlab("Cell type") +
  ylab("Expression") +
  ggtitle(genes[1]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(fill = F) +
  stat_compare_means(comparisons = my_compare)
ggsave("Fig2_HLA-E_mouse_CTC_barplot.pdf", width = 2.5, height = 4)
#-----




#' Figure 2
#' Figure 2I
#' barplot
#' KLRC1/KLRD1 expression levels in Blood
#-----
library(ggplot2)
library(ggpubr)

### 3 CD94 genes in blood
seu.immu <- readRDS(file = "../../P12345_immune_pc40_res2_plot2.rds")
seu.b <- subset(seu.immu, subset=group_tissue=="Portal-Venous")
table(seu.b$label)
Idents(seu.b) <- seu.b$label

genes <- c("KLRD1", "KLRC1")

# df to plot
df.plot <- as.data.frame(t(as.matrix(GetAssayData(seu.b, slot = "count"))[genes, ]))
head(df.plot)
df.plot$group <- Idents(seu.b)

## plot
ggbarplot(df.plot, x = "group", y=genes[1], add="mean_se",fill="group",  position=position_dodge(0.8),
          palette = pal_nejm()(6)) +
  xlab("Cell type") +
  ylab("Expression") +
  ggtitle(genes[1]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(fill = F)

ggsave("Fig2_CD94_blood_barplot.pdf", width = 4, height = 4)

ggbarplot(df.plot, x = "group", y=genes[2], add="mean_se",fill="group",  position=position_dodge(0.8),
          palette = pal_nejm()(6)) +
  xlab("Cell type") +
  ylab("Expression") +
  ggtitle(genes[2]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(fill = F)

ggsave("Fig2_NKG2A_blood_barplot.pdf", width = 4, height = 4)
#-----

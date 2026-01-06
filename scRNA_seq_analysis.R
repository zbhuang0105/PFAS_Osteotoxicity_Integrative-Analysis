library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)

# Note: Please place the 10X data (GSE147287 and GSE162454) in a "./data" subdirectory before running
data_dir <- paste0(getwd(),"/data")       
samples=list.files(data_dir)
dir=file.path(data_dir,samples)
afdata <- Read10X(data.dir = dir)
af <- CreateSeuratObject(counts = afdata, 
                         project = "SeuratObject", 
                         min.cells = 3,
                         min.features = 200)
afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = samples)
Idents(af)=afidens
af$Type=Idents(af)
af[["percent.mt"]] <- PercentageFeatureSet(af, pattern = "^MT-")
af[["percent.rb"]] <- PercentageFeatureSet(af, pattern = "^RP")
mask1 <- af$nCount_RNA >= 1000
mask2 <- af$nFeature_RNA >= 200 & af$nFeature_RNA <= 10000
mask3 <- af$percent.mt <= 20
mask4<-af$percent.rb<= 20
mask <-mask2 & mask3 & mask4
af <- af[, mask]



pdf(file = "01.vlnplot.pdf",width = 10,height = 7)
VlnPlot(af, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)+scale_fill_manual(values = c("#C77CFF","#7CAE00","#00BFC4","#F8766D"))
dev.off()


plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()

pdf(file = "01.corqc.pdf",width =9,height = 8)
plot2+plot3+plot_layout(ncol = 2)      
dev.off()


af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)


normalized.data <- af[["RNA"]]@data
normalized.data[1:20,1:4]
dim(normalized.data)


af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(af), 10)
top10


plot1 <- VariableFeaturePlot(af)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "01.topgene.pdf",width =7,height = 6)
plot2                   
dev.off()


af <- ScaleData(af)
scale.data <- af[["RNA"]]@scale.data
dim(scale.data)
scale.data[1:10,1:4]


all.genes <- rownames(af)
af <- ScaleData(af, features = all.genes)



af <- Seurat::RunPCA(af, features = VariableFeatures(object = af))
af <- Seurat::RunTSNE(af,dims = 1:20)
pdf(file = "02.rawtsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 3)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right") #top为图列位置最上方，除此之外还有right、left、bottom(意思同英文)
dev.off()
pdf(file = "02.rawpca.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "pca",pt.size = 3)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
colaa=distinctColorPalette(100)
pdf(file = "02.raw.tsne.split.pdf",width =12,height = 10)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =2
)
dev.off()


af <- RunHarmony(af, group.by.vars = "Type")

af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(af), 10)
top10


plot1 <- VariableFeaturePlot(af)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(file = "03.topgene.pdf",width =7,height = 6)
plot2                   
dev.off()



pdf(file = "03.harmony.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "harmony",pt.size = 3)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
af <- Seurat::RunTSNE(af,dims = 1:20,reduction ='harmony')
pdf(file = "03.tsne.pdf",width =7.5,height = 5.5)
DimPlot(af, reduction = "tsne",pt.size = 3)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
collist=c(ggsci::pal_nejm()(11))
names(collist)=names(table(af$Type))
pdf(file = "03.tsne.split.pdf",width =12,height = 7.5)
do_DimPlot(sample = af,
           plot.title = "",
           reduction = "tsne",
           legend.position = "bottom",
           dims = c(1,2),split.by = "Type",pt.size =2
) 
dev.off()


names(collist)=names(table(af$Type))

VizDimLoadings(af, dims = 1:2, reduction = "pca")

pdf(file = "03.pc_heatmap.pdf",width =7.5,height = 9)
DimHeatmap(af, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()

af <- JackStraw(af, num.replicate = 100)
af <- ScoreJackStraw(af, dims = 1:20)
pdf(file = "03.jackstrawplot.pdf",width =7.5,height = 5.5)
JackStrawPlot(af, dims = 1:20)
dev.off()
pdf(file = "03.ElbowPlot.pdf",width =5,height = 4)
ElbowPlot(af,ndims = 30)
dev.off()

afPC=19
af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
  af=FindClusters(af, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
apply(af@meta.data[,grep("RNA_snn_res",colnames(af@meta.data))],2,table)

p2_tree=clustree(af@meta.data, prefix = "RNA_snn_res.")
pdf(file = "03.clustertree.pdf",width =12,height =10)
p2_tree
dev.off()

af=FindNeighbors(af, dims = 1:afPC, reduction = "harmony")
af <- FindClusters(af, resolution = 1.5)


head(af@meta.data)
table(af@meta.data$seurat_clusters)
af.markers <- FindAllMarkers(af, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "04.cluster_markers.csv")
head(af.markers)


top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
col <- c(ggsci::pal_npg()(9),ggsci::pal_jco()(9),ggsci::pal_jama()(7),ggsci::pal_nejm()(8))
pdf(file = "04-cluster.hetmap.pdf",width =22,height = 16)
DoHeatmap(af,features = top5af.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()



af <- RunUMAP(af, dims = 1:afPC, reduction = "harmony")
af <- RunTSNE(af, dims = 1:afPC, reduction = "harmony")


pdf(file = "05-cluster.UMAP.pdf",width =6.5,height = 5.5)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "05-cluster.TSEN.pdf",width =6.5,height = 5.5)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


af <- FindClusters(af, resolution = 1.5)
genes <- list(
  "Osteoblastic OS cells" = c("COL1A1", "CDH11", "RUNX2"),
  "Proliferating osteoblastic OS cells" = c("COL1A1", "CDH11", "RUNX2", "TOP2A", "PCNA", "MKI67"),
  "Chondroblastic OS cells" = c("ACAN", "COL2A1", "SOX9"),
  "Osteoclastic cells" = c("CTSK", "MMP9"),
  "T/NK cells" = c("IL7R", "CD3D", "NKG7"),
  "Myeloid cells" = c("CD74", "CD14", "FCGR3A"),
  "Fibroblasts" = c("COL1A1", "LUM", "DCN"),
  "Pericytes" = c("ACTA2", "RGS5"),
  "MSCs" = c("CXCL12", "SFRP2", "MME"),
  "Myoblasts" = c("MYLPF", "MYL1"),
  "Endothelial cells" = c("PECAM1", "VWF")
)

pdf(file = "05.ann_cluster_marker.pdf",width =20,height = 30)
do_DotPlot(sample = af,features = genes,dot.scale = 12,colors.use = c("yellow","red"),legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

table(af@active.ident)
ann.ids <- c(
  "Myeloid cells",                     # cluster 0
  "Osteoblastic OS cells",            # cluster 1
  "Myeloid cells",            # cluster 2
  "Myeloid cells",                    # cluster 3
  "Osteoblastic OS cells",                    # cluster 4
  "Osteoblastic OS cells",            # cluster 5
  "Myeloid cells",              # cluster 6
  "Osteoblastic OS cells",            # cluster 7
  "Osteoclastic cells",            # cluster 8
  "Osteoclastic cells",                       # cluster 9
  "Osteoblastic OS cells",            # cluster 10
  "Fibroblasts",            # cluster 11
  "Osteoblastic OS cells",            # cluster 12
  "Osteoblastic OS cells",                    # cluster 13
  "MSCs",            # cluster 14
  "Osteoblastic OS cells",            # cluster 15
  "Osteoblastic OS cells",            # cluster 16
  "Myeloid cells",              # cluster 17
  "Myeloid cells",                    # cluster 18
  "Chondroblastic OS cells",               # cluster 19
  "Endothelial cells",         # cluster 20
  "Myeloid cells",                    # cluster 21
  "Fibroblasts",                       # cluster 22
  "Osteoblastic OS cells",                    # cluster 23
  "Osteoblastic OS cells",                        # cluster 24
  "Osteoblastic OS cells",                      # cluster 25
  "Myoblasts",            # cluster 26
  "Osteoclastic cells",                    # cluster 27
  "Fibroblasts",            # cluster 28
  "Osteoblastic OS cells",                      # cluster 29
  "Osteoblastic OS cells",              # cluster 30
  "Osteoblastic OS cells",            # cluster 31
  "Osteoclastic cells",              # cluster 32
  "Osteoclastic cells", # cluster 33
  "Myeloid cells",            # cluster 34
  "Myeloid cells",                    # cluster 35
  "Osteoclastic cells",              # cluster 36
  "Osteoclastic cells",                    # cluster 37
  "Myeloid cells",            # cluster 38
  "Osteoblastic OS cells"                     # cluster 39
)


afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens
af$cellType=Idents(af)

# 可视化UMAP/tSNE
pdf(file = "05-ann.scRNA.UMAP.pdf",width =15,height = 10)
DimPlot(af, reduction = "umap", label = T, label.size = 3.5,pt.size = 2.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
pdf(file = "05-ann.scRNA.TSEN.pdf",width =10,height = 10)
DimPlot(af, reduction = "tsne", label = T, label.size = 3.5,pt.size = 2.5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()


pdf(file = "05-ann.scRNA.UMAP.pdf", width = 20, height = 15)
DimPlot(
  af, reduction = "umap", label = TRUE, label.size = 8.0, pt.size = 6.0
) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    legend.position = "right",
    axis.text = element_text(size = 20),      
    axis.title = element_text(size = 22),     
    legend.text = element_text(size = 18),    
    legend.title = element_text(size = 20)    
  )
dev.off()


pdf(file = "05-ann.scRNA.TSNE.pdf", width = 17, height = 15)
DimPlot(
  af, reduction = "tsne", label = TRUE, label.size = 8.0, pt.size = 6.0
) +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )
dev.off()

af.markers <- FindAllMarkers(af, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(af.markers,file = "06.cell_markers.csv")
# get top 10 genes
top5af.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# plot
pdf(file = "06-cell_marker.hetmap.pdf",width =10,height = 10)
DoHeatmap(af,features = top5af.markers$gene,
          group.colors = col) +
  ggsci::scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')
dev.off()

df=as.data.frame(af.markers[,c("gene","cluster")])
table(df$cluster)
dirnb=(paste0(getwd(),"/D1.afgmt.gmt"))
gmtlist=list()
for (afk in as.character(names(table(df$cluster)))) {
  ddgene=rownames(df)[df$cluster==afk]
  aaname=afk
  gmtlist[[aaname]]=ddgene
}
output_gmt<-function(geneset,file){
  sink(file) 
  lapply(names(geneset),function(i){
    cat(paste(c(i,'NA',geneset[[i]]),collapse='\t')) 
    cat('\n') 
  })
  sink() 
}
output_gmt(gmtlist,dirnb)


colaa=distinctColorPalette(100)

afgenes=c("AKT1",    
          "ESR1",
          "HSP90AA1",
          "PARP1",
          "EGFR",
          "CASP3",
          "ALB",
          "ANXA5",
          "SRC",
          "MDM2")  
pdf(file = "06-cell_FeaturePlot.pdf",width =12,height = 6)
FeaturePlot(af, features = afgenes, cols = c("grey", "red"),min.cutoff = 1, max.cutoff = 3,ncol=4,pt.size = 0.5)    #min.cutoff与max.cutoff修改截断以更好可视化结果，通过颜色强调基因的分布
dev.off()
pdf(file = "06-cell_VlnPlot.pdf",width =8,height = 5)
VlnPlot(af, features = afgenes,group.by = "cellType", stack=TRUE,cols = colaa)+ NoLegend()   
dev.off()

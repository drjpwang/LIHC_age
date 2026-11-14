rm(list = ls())
setwd("/data/LIHC/")
if(! dir.exists("./16_Hepatocytes")){
  dir.create("./16_Hepatocytes")
}
setwd("./16_Hepatocytes")

getwd()

library(Seurat)

load("../05_scRNA-seq/LIHC_sc_anno.rda")
UMAP <- sc_dat2
head(UMAP@meta.data)
table(UMAP@meta.data$cell_type)
dim(UMAP)

sc_anno <- DimPlot(UMAP, reduction = "umap", label = TRUE, pt.size = 3, split.by = "type")


cnvcluster <- UMAP[,Idents(UMAP) %in% c("Hepatocytes","Fibroblast")]
dim(cnvcluster)
# [1] 32695 18899

#数据拆分
head(cnvcluster@meta.data, 5)
table(cnvcluster$orig.ident)
cnvcluster[["RNA"]] <- split(cnvcluster[["RNA"]], f = cnvcluster$orig.ident)

#标准化,将每个基因的counts除以细胞总的counts数，乘上10,000，再进行对数转换
cnvcluster <- NormalizeData(cnvcluster)
#特征选择,选择在细胞与细胞之间具有高度变异性的基因,默认选择前2000个高度变异基因
cnvcluster <- FindVariableFeatures(cnvcluster, selection.method = "vst", nfeatures = 2000)

cnvcluster <- ScaleData(cnvcluster)
# cnvcluster <- RunPCA(cnvcluster)

#使用harmony方法去批次
# cnvcluster1 <- IntegrateLayers(
#   object = cnvcluster, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "harmony",
#   verbose = FALSE
# )

cnvcluster1 <- JoinLayers(cnvcluster)

#线性PCA降维
# aa_sc1 <- RunPCA(aa_sc1, features = VariableFeatures(object = aa_sc1))
# sc_pca1 <- DimPlot(cnvcluster1, reduction = "cca")
# ggsave("01_PCA_plot.pdf", sc_pca1, width = 8, height = 6)
# ggsave("01_PCA_plot.png", sc_pca1, width = 8, height = 6)

#DimHeatmap(luad_data1, dims = 1:5, cells = 500, balanced = TRUE)

#维数选择
# cnvcluster1 <- JackStraw(cnvcluster1, num.replicate = 100, dims = 20)
# cnvcluster1 <- ScoreJackStraw(cnvcluster1, dims = 1:20)

# sc_jack <- JackStrawPlot(cnvcluster1, dims = 1:20)
# ggsave("02_JackStraw_plot.pdf", sc_jack, width = 8, height = 6)
# ggsave("02_JackStraw_plot.png", sc_jack, width = 8, height = 6)

# sc_elbow <- ElbowPlot(cnvcluster1, ndims = 20)
# ggsave("03_Elbow_plot.pdf", sc_elbow, width = 14, height = 8)
# ggsave("03_Elbow_plot.png", sc_elbow, width = 14, height = 8)

#使用clustree测试cluster分类效果确定最合适的分辨率resolution
# library(clustree)
# cnvcluster2 <- cnvcluster1
# cnvcluster2 <- FindNeighbors(cnvcluster2,dim=1:10,reduction = "cca")
# cnvcluster2 <- FindClusters(cnvcluster2, resolution = c(seq(0.1,0.9,.1)))
#length(levels(Idents(luad_data2)))
# clustree(cnvcluster2@meta.data, prefix = "RNA_snn_res.")

# ggsave("clustree.pdf", height = 10, width = 6)

#细胞聚类
# dims = 1:10
# cnvcluster1 <- FindNeighbors(cnvcluster1, dims = dims, reduction = "cca")
# cnvcluster1 <- FindClusters(cnvcluster1, resolution = 0.1)
# table(Idents(cnvcluster1))

#非线性降维（UMAP/tSNE）
# cnvcluster1 <- RunUMAP(cnvcluster1, dims = dims, reduction = "cca")
# sc_raw <- DimPlot(cnvcluster1, reduction = "umap", label = TRUE, pt.size = 5)

# ggsave("04_raw_umap.pdf", sc_raw, width = 10, height = 8)
# ggsave("04_raw_umap.png", sc_raw, width = 10, height = 8)

# cnvcluster1 <- RunTSNE(cnvcluster1, dims = dims, reduction = "cca")
# sc_raw1 <- DimPlot(macro_dat1, reduction = "tsne", label = TRUE, pt.size = 5)

# ggsave("05_raw_tsne.pdf", sc_raw1, width = 10, height = 8)
# ggsave("05_raw_tsne.png", sc_raw1, width = 10, height = 8)


# marker <- c("EPCAM","KRT19","KRT18","CDH1", # epithelial"
#             "DCN","THY1","COL1A1","COL1A2", # Fibroblast
#             "PECAM1","VWF","CLDN5","RAMP2") # endothelial

# sc_dot <- DotPlot(cnvcluster1, features = marker)+
#   coord_flip() + #翻转
#   theme(panel.grid = element_blank(), 
#         axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
#   labs(x=NULL,y=NULL) + 
#   guides(size = guide_legend("Percent Expression") )+ #legend
#   scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
# sc_dot
# ggsave("05_dotplot_raw.pdf", sc_dot, width = 10, height = 8)
# ggsave("05_dotplot_raw.png", sc_dot, width = 10, height = 8)

# macro_marker <- FindAllMarkers(macro_dat1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# macro_Htop20 <- macro_marker %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
# ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少
# sc_heat <- DoHeatmap(macro_dat1, features = macro_Htop20$gene) + NoLegend()
# sc_heat
# # 
# ggsave("06_heatmap.pdf", sc_heat, width = 10, height = 8)
# ggsave("06_heatmap.png", sc_heat, width = 10, height = 8)

# save(cnvcluster1, file = "scRNA_anno_before.rda")

# sc_type <- DimPlot(macro_dat1, reduction = "umap", group.by = "orig.ident", pt.size = 2)
# ggsave("07_batch_umap_raw.pdf", sc_type, width = 10, height = 8)
# ggsave("07_batch_umap_raw.png", sc_type, width = 10, height = 8)

#细胞注释
# sc_dat2 <- cnvcluster
# sc_dat2 <- RenameIdents(sc_dat2, '0'='Epithelial cells', '1'='Endothelial cells', '2'='Fibroblasts',
#                         '3'='Fibroblasts', '4'='Epithelial cells', '5'='Epithelial cells',
#                         '6'='Epithelial cells', '7'='Epithelial cells')
# sc_dat2@meta.data$cell_type <- Idents(sc_dat2)
# sc_anno <- DimPlot(sc_dat2, reduction = "umap", label = TRUE, pt.size = 5) +
#   theme(text=element_text(size= 13 ))#+NoLegend()
# ggsave("06_anno_umap.pdf", sc_anno, width = 14, height = 10)
# ggsave("06_anno_umap.png", sc_anno, width = 14, height = 10)

# sc_anno1 <- DimPlot(sc_dat2, reduction = "umap", label = TRUE, pt.size = 5, split.by = "Group") +
#   theme(text=element_text(size= 13 ))#+NoLegend()
# sc_anno1


# sc_dot1 <- DotPlot(sc_dat2, features = marker)+
#   coord_flip() + #翻转
#   theme(panel.grid = element_blank(), 
#         axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
#   labs(x=NULL,y=NULL) + 
#   guides(size = guide_legend("Percent Expression") )+ #legend
#   scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
# sc_dot1
# ggsave("07_dotplot_anno.pdf", sc_dot1, width = 10, height = 8)
# ggsave("07_dotplot_anno.png", sc_dot1, width = 10, height = 8)
# sc_type <- DimPlot(luad_data3, reduction = "umap", group.by = "project")

# save(sc_dat2, file = "scRNA_sc_anno.rda")

# head(sc_dat2@meta.data)
# table(sc_dat2$Group)


# sc_dat2 <- cnvcluster1
#识别肿瘤细胞
#使用inferCNV区分肿瘤细胞和正常细胞,需要三个文件
#1.基因表达矩阵（未标准化的原count矩阵），每一行为基因，每一列为细胞(Raw Counts Matrix)
# dfcount = as.data.frame(sc_dat2@assays$RNA$counts)
# dfcount <- cnvcluster@assays$RNA$count
dfcount <- GetAssayData(cnvcluster1, slot = "count", assay = "RNA")
dfcount <-as.data.frame(dfcount)
dfcount[1:3,1:3]
dfcount1 <- dfcount


#2.细胞注释文件包括两列，第一列为细胞名称，第二列为其细胞类型
groupinfo= data.frame(cellID = colnames(dfcount))
groupinfo$cellType = as.vector(cnvcluster$cell_type)
groupinfo[1:6,]
#给参考的非肿瘤细胞提供标签(non-malignant)
groupinfo1 <- groupinfo
table(groupinfo1$cellType)
library(stringr)
for (i in 1:dim(groupinfo1)[1]) {
  if(groupinfo1[i,2]=="Fibroblast"){
    groupinfo1[i,2]=str_c(groupinfo1[i,2],"(non-malignant)",sep = '')
  }
  # else
  #   if(groupinfo1[i,2]=="Fibroblasts"){
  #     groupinfo1[i,2]=str_c(groupinfo1[i,2],"(non-malignant)",sep = '')
  #   }
}
groupinfo1[1:7,]
table(groupinfo1$cellType)
#测试
library(stringr)
# groupinfo[5,2] <- str_c(groupinfo[5,2],"(non-malignant)",sep = '')


#3.基因坐标文件，是每个基因在基因组上位置的注释
genenames <- rownames(dfcount)

#for (i in 1:length(genenames)) {

#genenames[i] <- strsplit(genenames[i], "\\.")[[1]][1]

#}
#test <- genenames[1:5]
#for (i in 1:5) {

#   test[i] <- strsplit(test[i], "\\.")[[1]][1]
#}
library(AnnoProbe)
geneInfor = annoGene(genenames,"SYMBOL",'human')
head(geneInfor)
geneInfor1 <- geneInfor
geneInfor1 <- geneInfor1[,c(1,4,5,6)]
#setdiff可以取两者中不同的元素，只取后者不在前者中的
length(setdiff(genenames, geneInfor1[,1]))
geneInfor1 <- geneInfor1[!duplicated(geneInfor1[,1]),]

#取出成功注释的gene
#dfcount1$genenames <- rownames(dfcount1)
#dfcount1 <- dfcount1[!duplicated(dfcount1$genenames),]
dfcount1 <- dfcount1[geneInfor1$SYMBOL,]

#储存为文档
write.table(dfcount1, "expfile.txt",sep = "\t", quote = F)
write.table(groupinfo1, "groupfile.txt", sep = "\t",quote=F,row.names = F,col.names = F)
write.table(geneInfor1, "genefile.txt",sep='\t',quote=F,row.names = F,col.names = F)

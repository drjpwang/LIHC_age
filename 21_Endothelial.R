rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./21_Endothelial")){
  dir.create("./21_Endothelial")
}
setwd("./21_Endothelial")

getwd()


#筛选关键差异细胞
load("../19_scRNA-seq_inferCNV/LIHC_sc_anno.rda")
table(sc_dat2$cell_type)

#选取内皮细胞作为关键细胞
endo_dat <- sc_dat2[,Idents(sc_dat2) %in% "Endothelial"]
dim(endo_dat)
# [1] 32695  4478
#数据拆分
head(endo_dat@meta.data, 5)
endo_dat[["RNA"]] <- split(endo_dat[["RNA"]], f = endo_dat$orig.ident)

#标准化,将每个基因的counts除以细胞总的counts数，乘上10,000，再进行对数转换
endo_dat <- NormalizeData(endo_dat)
#特征选择,选择在细胞与细胞之间具有高度变异性的基因,默认选择前2000个高度变异基因
endo_dat <- FindVariableFeatures(endo_dat, selection.method = "vst", nfeatures = 2000)

endo_dat <- ScaleData(endo_dat)
endo_dat <- RunPCA(endo_dat)

endo_dat1 <- endo_dat
#使用cca方法去批次
# endo_dat1 <- IntegrateLayers(
#   object = endo_dat, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "cca",
#   verbose = FALSE
# )

endo_dat1 <- JoinLayers(endo_dat1)

#线性PCA降维
# aa_sc1 <- RunPCA(aa_sc1, features = VariableFeatures(object = aa_sc1))
sc_pca1 <- DimPlot(endo_dat1, reduction = "pca")
ggsave("01_PCA_plot.pdf", sc_pca1, width = 8, height = 6)
ggsave("01_PCA_plot.png", sc_pca1, width = 8, height = 6)

#DimHeatmap(luad_data1, dims = 1:5, cells = 500, balanced = TRUE)

#维数选择
endo_dat1 <- JackStraw(endo_dat1, num.replicate = 100)
endo_dat1 <- ScoreJackStraw(endo_dat1, dims = 1:20)

sc_jack <- JackStrawPlot(endo_dat1, dims = 1:20)
ggsave("02_JackStraw_plot.pdf", sc_jack, width = 8, height = 6)
ggsave("02_JackStraw_plot.png", sc_jack, width = 8, height = 6)

sc_elbow <- ElbowPlot(endo_dat1)
ggsave("03_Elbow_plot.pdf", sc_elbow, width = 14, height = 8)
ggsave("03_Elbow_plot.png", sc_elbow, width = 14, height = 8)

#使用clustree测试cluster分类效果确定最合适的分辨率resolution
library(clustree)
endo_dat2 <- endo_dat1
endo_dat2 <- FindNeighbors(endo_dat2,dim=1:15,reduction = "harmony")
endo_dat2 <- FindClusters(endo_dat2, resolution = c(seq(0.1,0.9,.1)))
#length(levels(Idents(luad_data2)))
clustree(endo_dat2@meta.data, prefix = "RNA_snn_res.")

#细胞聚类
dims = 1:15
endo_dat1 <- FindNeighbors(endo_dat1, dims = dims, reduction = "harmony")
endo_dat1 <- FindClusters(endo_dat1, resolution = 0.2)
table(Idents(endo_dat1))

#非线性降维（UMAP/tSNE）
endo_dat1 <- RunUMAP(endo_dat1, dims = dims, reduction = "harmony")
sc_raw <- DimPlot(endo_dat1, reduction = "umap", label = TRUE, pt.size = 5)

ggsave("04_raw_umap.pdf", sc_raw, width = 10, height = 8)
ggsave("04_raw_umap.png", sc_raw, width = 10, height = 8)

endo_dat1 <- RunTSNE(endo_dat1, dims = dims, reduction = "harmony")
sc_raw1 <- DimPlot(endo_dat1, reduction = "tsne", label = TRUE, pt.size = 5)

ggsave("05_raw_tsne.pdf", sc_raw1, width = 10, height = 8)
ggsave("05_raw_tsne.png", sc_raw1, width = 10, height = 8)

# mono_marker <- FindAllMarkers(endo_dat1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# mono_Htop10 <- mono_marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少
# sc_heat <- DoHeatmap(mono_dat1, features = mono_Htop10$gene) + NoLegend()
# 
# ggsave("06_heatmap.pdf", sc_heat, width = 10, height = 8)
# ggsave("06_heatmap.png", sc_heat, width = 10, height = 8)

save(endo_dat1, file = "Endotheial_anno_before.rda")

#细胞注释
#从文献中获取marker基因
marker <- c("CCL21","MMRN1","PROX1", # lymphatic ECs 淋巴内皮
            "GLUL","FBLN5","HEY1","IGFBP3", # Arteries 动脉
            "GJA4","BGN",
            "ACKR1","AHNAK","GPIHBP1","CAVIN2",          # Veins 静脉
            "OIT3","CLEC4G","CLEC1B",   # "CD34","CDH5",
            "IL6ST","SGK1","CD36","CD4")                     # Capillaries 毛细血管


DotPlot(endo_dat1, features = marker)
sc_dot <- DotPlot(endo_dat1, features = marker, dot.scale = 10)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
sc_dot
ggsave("07_dotplot_raw.pdf", sc_dot, width = 10, height = 8)
ggsave("07_dotplot_raw.png", sc_dot, width = 10, height = 8)

sc_type <- DimPlot(endo_dat1, reduction = "umap", group.by = "orig.ident", pt.size = 2)
ggsave("08_batch_umap_raw.pdf", sc_type, width = 10, height = 8)
ggsave("08_batch_umap_raw.png", sc_type, width = 10, height = 8)

#细胞注释
sc_dat2 <- endo_dat1
sc_dat2 <- RenameIdents(sc_dat2, '0'='Vein ECs', '1'='Capillary ECs', 
                        '2'='Lymphatic ECs', '3'='Artery-1 ECs',
                        '4'='Artery-2 ECs', '5'='Capillary ECs')
sc_dat2@meta.data$cell_type <- Idents(sc_dat2)
sc_anno <- DimPlot(sc_dat2, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave("09_anno_umap.pdf", sc_anno, width = 10, height = 8)
ggsave("09_anno_umap.png", sc_anno, width = 10, height = 8)

sc_dot1 <- DotPlot(sc_dat2, features = marker, dot.scale = 10)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
ggsave("10_dotplot_anno.pdf", sc_dot1, width = 10, height = 8)
ggsave("10_dotplot_anno.png", sc_dot1, width = 10, height = 8)
# sc_type <- DimPlot(luad_data3, reduction = "umap", group.by = "project")

save(sc_dat2, file = "Endothelial_sc_anno.rda")


#计算差异基因
endo.markers <- FindAllMarkers(sc_dat2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(endo.markers, "Endothelial_FindAllMarkers.csv")
endo_marker <- unique(endo.markers$gene)

#读取hubgene
load("../03_Survival/LASSO_GeneID.rda")

endo_marker <- union(endo_marker, lasso_geneids$symbol)
setdiff(lasso_geneids$symbol,endo.markers$gene)
save(endo_marker, file = "endo_markers.rda")
write.csv(endo_marker, "endo_markers.csv")



if(! dir.exists("./Gene_expression")){
  dir.create("./Gene_expression")
}
setwd("./Gene_expression")

getwd()

load("../03_Survival/LASSO_GeneID.rda")

feaplot <- FeaturePlot(sc_dat2,
                       features = lasso_geneids$symbol,repel = T, label = F,
                       slot = "data", cols = c("lightgrey","red"), ncol = 2, order = T)
ggsave("03_feature_featureplot.pdf", feaplot, width = 10, height = 10)
ggsave("03_feature_featureplot.png", feaplot, width = 10, height = 10)

vlplot <- VlnPlot(sc_dat2, features = lasso_geneids$symbol, group.by = 'cell_type', pt.size = 0)
ggsave("01_feature_vlnplot.pdf", vlplot, width = 10, height = 8)
ggsave("01_feature_vlnplot.png", vlplot, width = 10, height = 8)

dotplot <- DotPlot(sc_dat2, features = lasso_geneids$symbol)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
ggsave("02_feature_dotplot.pdf", dotplot, width = 10, height = 8)
ggsave("02_feature_dotplot.png", dotplot, width = 10, height = 8)


scRNA.norm.umap <- sc_dat2
table(scRNA.norm.umap$cell_type)
sce.all =  AddModuleScore(object = scRNA.norm.umap,list(lasso_geneids$symbol))
head(sce.all@meta.data)
stat.test <- sce.all@meta.data%>%
  group_by(cell_type)%>%
  wilcox_test(Cluster1 ~ type)%>%
  adjust_pvalue(method = 'fdr')

stat.test$p <- ifelse(stat.test$p<0.001,"***",
                      ifelse(stat.test$p<0.01,"**",
                             ifelse(stat.test$p<0.05,"*",'ns')))
temp1 <- stat.test[stat.test$p < 0.05, ]
write.csv(temp1,"temp1.csv")
table(sce.all@meta.data$sample_type)

External_va <- ggplot(sce.all@meta.data, aes(x = cell_type, y = Cluster1, color = type)) +
  # ylim(0,0.8)+
  stat_boxplot(geom="errorbar",
               width=0.2,
               position = position_dodge(0.9)) +
  geom_boxplot(aes(color=type,fill=type),
               width=0.8,
               position=position_dodge(0.9),
               outlier.shape = NA  )+
  scale_fill_manual(values= c("#4682B4","#CD3700"))+#绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  #
  # annotate("rect", xmin = 0.4, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#c1f1fc") +
  #  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#ebffac") +
  #添加散点
  # geom_dotplot(dotsize = 0.4,binaxis = "y", stackdir = "center",position = position_dodge(0.8))+
  #添加分组辅助线
  #geom_vline(xintercept = 1.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  #geom_vline(xintercept = 2.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  #geom_vline(xintercept = 3.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  #geom_vline(xintercept = 4.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12),
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title=" ", x="", y = "AddModuleScore",size=20) +
  scale_colour_manual(values = c("black","black"))+
  #   scale_color_brewer(palette = 'Set2')+
  stat_compare_means(aes(group = type,label = "p.format"),
                     method = "wilcox.test",label = 'p.signif')#+,label.x = 1.4
# stat_pvalue_manual(stat.test,
#             y.position = 0.2,
#            size = 3.2,
#            family = "Times",
#            label = "p",
#            parse = T,
#            face = "bold")

External_va
ggsave(filename = '04_AddModuleScore_diff.pdf', External_va, w=8, h=6)
ggsave(filename = '04_AddModuleScore_diff.png', External_va, w=8, h=6, dpi = 600)
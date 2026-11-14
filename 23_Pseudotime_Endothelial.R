rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./23_Pseudotime_Endothelial")){
  dir.create("./23_Pseudotime_Endothelial")
}
setwd("./23_Pseudotime_Endothelial")

getwd()


#拟时序分析
library(monocle)

load("../21_Endothelial/Endothelial_sc_anno.rda")
#使用count
epicount1 <- GetAssayData(sc_dat2, slot = "count", assay = "RNA")
#细胞对象，表型信息
pd <- new('AnnotatedDataFrame', data = sc_dat2@meta.data)
#基因对象
fdata <- data.frame(gene_short_name = row.names(epicount1), row.names = row.names(epicount1))
fd <- new('AnnotatedDataFrame', data = fdata)

#构建S4对象，CellDataSet
#单细胞稀疏矩阵的话用negbinomial.size()，如果是UMI的话不要标准化；
#FPKM值用tobit()；logFPKM值用gaussianff()
est <- newCellDataSet(epicount1, 
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

#估计size factors和 dispersions
est <- estimateSizeFactors(est)
est <- estimateDispersions(est)
#过滤低质量细胞,过滤低于1%细胞中检出的基因
est <- detectGenes(est, min_expr = 0.1)#计算每个基因在多少细胞中表达

head(fData(est))
#expressed_genes <- row.names(subset(fData(est),num_cells_expressed >= 10))
head(pData(est))
dim(est)

#选择进行拟时序的基因
#使用FindAllmarkers之后的各个亚群的marker基因
load("../21_Endothelial/endo_markers.rda")
pc_genes <- unique(endo_marker)
# epi_genes <- epi_genes[1:1500]

#读取关键基因
load("../03_Survival/LASSO_GeneID.rda")
sur_g <- lasso_geneids$symbol

pc_genes <- union(pc_genes, sur_g)
setdiff(lasso_geneids$symbol, pc_genes)

test_ordering_genes <- pc_genes
est <- setOrderingFilter(est, ordering_genes = test_ordering_genes)

plot_ordering_genes(est)

#推断轨迹，并按照拟时序给细胞排序
#降维和排序
#减少其他因素的影响，比如不同样本、不同批次
est <- reduceDimension(est, max_components = 2, reduction_method = 'DDRTree', 
                       residualModelFormulaStr = "~num_genes_expressed")
est <- orderCells(est)#, root_state = 5

#报以下错时输入以下语句进行源代码修改
#Error in if (class(projection) != "matrix") projection <- as.matrix(projection) : 
#the condition has length > 1
#trace('project2MST', edit = T, where = asNamespace("monocle"))

save(est, file = "Pseudotime_result.rda")


#绘图
plot1 <- plot_cell_trajectory(est,color_by = "Pseudotime")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 16))
# ggsave("01_monocle_pseudotime.pdf", plot1, width = 10, height = 10)
plot2 <- plot_cell_trajectory(est,color_by = "cell_type")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 8))
# ggsave("02_monocle_seurat_clusters.pdf", plot2, width = 14, height = 8)
plot3 <- plot_cell_trajectory(est,color_by = "State")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
# ggsave("03_monocle_state.pdf", plot3, width = 14, height = 8)
plot4 <- plot_cell_trajectory(est,color_by = "type")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

plot1 + plot3
ggsave("01_monocle_pseudotime&state.pdf", plot1+plot3, width = 18, height = 8)
ggsave("01_monocle_pseudotime&state.png", plot1+plot3, width = 25, height = 8)

plot2
ggsave("01_monocle_cell_type.pdf", plot2, width = 18, height = 8)
ggsave("01_monocle_cell_type.png", plot2, width = 25, height = 8)

plot4
ggsave("01_monocle_sample_type.pdf", plot4, width = 12, height = 8)
ggsave("01_monocle_sample_type.png", plot4, width = 12, height = 8)

# plot2 + plot3
# ggsave("05_monocle_seurat_clusters&state.pdf", plot2+plot3, width = 18, height = 8)
# 
# plot1 + plot3
# ggsave("06_monocle_pseudotime&state.pdf", plot1+plot3, width = 18, height = 8)
# 
# #树形图
# plot5 <- plot_complex_cell_trajectory(est, x=1, y=2, color_by = "seurat_clusters")
# ggsave("07_cell_trajectory_seurat_clusters.pdf", plot5, width = 18, height = 8)
# 
# plot6 <- plot_cell_trajectory(est, color_by = "seurat_clusters") +
#   facet_wrap(~seurat_clusters, nrow = 2)
# ggsave("08_cell_trajectory_seurat_clusters_wrap.pdf", plot6, width = 18, height = 8)

# saveRDS(est, 'CA_pseudotime_result.RDS')

#est <- readRDS('macro_pseudotime_result.RDS')
#绘制画沿时间轴的细胞密度图
library(ggpubr)
library(ggplot2)
#pData(est)取出的是est对象中est@phenoData的内容
df <- pData(est)
plot7 <- ggplot(df, aes(Pseudotime, color = cell_type, fill = cell_type))+
  geom_density(bw=0.5, size=1, alpha=0.5) + theme_classic2()
ggsave("02_monocle_density.pdf", plot7, width = 10, height = 8)
ggsave("02_monocle_density.png", plot7, width = 10, height = 8)

# #自己设置颜色进行展示
# density_col <- c("Vascular smooth muscle contraction"="#DC143C","Cell cycle"="#0000FF","Immune"="#00FF4DFF",
#                  "Cell adhesion molecules"="#CD00CD","Ribosome"="#FF1493","Cell cycle&Immune"="#9370DB",
#                  "Basal"="#A6D854","Luminal1"="#FF6347","Luminal2"="#FFFF00","Undefined"="#87CEFA")
# plot9 <- ggplot(df, aes(Pseudotime, color = seurat_clusters, fill = seurat_clusters))+
#   geom_density(bw=0.5, size=1, alpha=0.5) + theme_classic2()+
#   scale_fill_manual(name="",values=density_col)+scale_color_manual(name="",values = density_col)

#查看生物标志物拟时序表达情况
biomarker <- sur_g

plot8 <- plot_genes_in_pseudotime(est[biomarker,],
                                  color_by = "cell_type",
                                  ncol = 2)+
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
ggsave("03_Biomarker.pdf", plot8, width = 10, height = 8)
ggsave("03_Biomarker.png", plot8, width = 10, height = 8)


trace("plot_pseudotime_heatmap", edit = T)
p_pseudotime <- plot_pseudotime_heatmap(est[biomarker,], 
                                        cores = 6, 
                                        num_clusters = 6, 
                                        show_rownames = T, return_heatmap = T,
                                        hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
ggsave("04_Biomarker_heatmap.pdf", p_pseudotime, width = 10, height = 8)
ggsave("04_Biomarker_heatmap.png", p_pseudotime, width = 10, height = 8)

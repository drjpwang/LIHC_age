rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./05_scRNA-seq")){
  dir.create("./05_scRNA-seq")
}
setwd("./05_scRNA-seq")

getwd()


library(Seurat)

sa_sc1 <- readRDS("~/LIHC/analysis_aging/01.Seurat/merge/harmony_0.3_20/Seurat.rds")

#使用clustree测试cluster分类效果确定最合适的分辨率resolution
library(clustree)
sa_sc2 <- sa_sc1
sa_sc2 <- FindNeighbors(sa_sc2,dim=1:30,reduction = "harmony")
sa_sc2 <- FindClusters(sa_sc2, resolution = c(seq(0.1,0.9,.1)))
#length(levels(Idents(luad_data2)))
library(clustree)
clustree(sa_sc2@meta.data, prefix = "RNA_snn_res.")

ggsave("clustree_plot.pdf", width = 10, height = 8)

#细胞聚类
dims = 1:30
sa_sc1 <- FindNeighbors(sa_sc1, dims = dims, reduction = "harmony")
sa_sc1 <- FindClusters(sa_sc1, resolution = 0.6)
table(Idents(sa_sc1))

#非线性降维（UMAP/tSNE）
sa_sc1 <- RunUMAP(sa_sc1, dims = dims, reduction = "harmony")
sc_raw <- DimPlot(sa_sc1, reduction = "umap", label = TRUE, pt.size = 3)

ggsave("06_raw_umap.pdf", sc_raw, width = 10, height = 8)
ggsave("06_raw_umap.png", sc_raw, width = 10, height = 8)

sa_sc1 <- RunTSNE(sa_sc1, dims = dims, reduction = "harmony")
sc_raw1 <- DimPlot(sa_sc1, reduction = "tsne", label = TRUE, pt.size = 3)

ggsave("07_raw_tsne.pdf", sc_raw1, width = 10, height = 8)
ggsave("07_raw_tsne.png", sc_raw1, width = 10, height = 8)

#细胞注释
#从文献中获取marker基因
marker <- c("EPCAM","KRT19","CDH1", # epithelial "KRT18",
            "DCN","THY1","COL1A1","COL1A2", # fibroblast
            "PECAM1","VWF","CLDN5","RAMP2", # endothelial
            # "ACTA2","TAGLN",         # Smooth Muscle cell
            "ALB","APOA2","TTR",            # Hepatocytes
            # "RGS5","MYL9","KCNJ8",           # Pericytes
            # "ITGAX",                         # Myeloid
            "CD3D","CD3E","CD3G","TRAC",    # T cell
            "NKG7","GNLY","KLRD1",        # NK cell
            "IGHM","IGHG3","IGHA2","CD79A", # B cell
            "MS4A1","CD19",      # Plasma
            "LYZ","CLEC9A","FCER1A",        # DC  "CCR7",
            # "CSF3R","FCGR3B",               # Neutrophil
            # "KIT","MS4A2","GATA2",          # Mast cell
            "CD163","CCL18","CD68")  # macrophage  ,"CD86"
# cell <- c("epithelial","epithelial","epithelial","Malignant cell","Malignant cell","Malignant cell",
#           "fibroblast","fibroblast","fibroblast","fibroblast",
#           "endothelial","endothelial","endothelial","endothelial",
#           "Smooth Muscle cell","Smooth Muscle cell","Pericytes","Pericytes","Pericytes",
#           "Myeloid","Myeloid","T cell","T cell","T cell","T cell",
#           "B cell","B cell","B cell","B cell","Plasma","NK cell","NK cell","NK cell")
# marker_gene <- data.frame(marker=marker, celltype=cell)
# write.csv(marker_gene, "cell_markergene.csv")
sc_dot <- DotPlot(sa_sc1, features = marker)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
sc_dot
ggsave("08_dotplot_raw.pdf", sc_dot, width = 10, height = 8)
ggsave("08_dotplot_raw.png", sc_dot, width = 10, height = 8)

sc_type <- DimPlot(sa_sc1, reduction = "umap", group.by = "orig.ident")
ggsave("09_batch_umap_raw.pdf", sc_type, width = 10, height = 8)
ggsave("09_batch_umap_raw.png", sc_type, width = 10, height = 8)

save(sa_sc1, file = "LIHC_sc_before_anno.rda")

#细胞注释
sc_dat2 <- sa_sc1
sc_dat2 <- RenameIdents(sc_dat2,
                        '0'='T cells', '1'='NKT cells', '2'='T cells',  
                        '3'='T cells', '4'='Macrophage', '5'='NK cells',
                        '6'='Endothelial', '7'='Hepatocytes', '8'='T cells',
                        '9'='T cells', '10'='Macrophage', '11'='Plasma cells',
                        '12'='Macrophage', '13'='Fibroblast', '14'='B cells',
                        '15'='Hepatocytes', '16'='T cells', '17'='Epithelial', 
                        '18'='DC','19'='Macrophage', '20'='Hepatocytes', 
                        '21'='Hepatocytes', '22'='Macrophage', '23'='Hepatocytes',
                        '24'='Hepatocytes', '25'='NK cells', '26'='Hepatocytes',
                        '27'='Hepatocytes', '28'='Macrophage', '29'='Hepatocytes',
                        '30'='Macrophage', '31'='T cells', '32'='Macrophage',
                        '33'='Hepatocytes', '34'='Hepatocytes', '35'='NKT cells',
                        '36'='T cells')
sc_dat2@meta.data$cell_type <- Idents(sc_dat2)
sc_anno <- DimPlot(sc_dat2, reduction = "umap", label = TRUE, pt.size = 1)
sc_anno
ggsave("10_anno_umap.pdf", sc_anno, width = 10, height = 8)
ggsave("10_anno_umap.png", sc_anno, width = 10, height = 8)

sc_anno1 <- DimPlot(sc_dat2, reduction = "tsne", label = TRUE, pt.size = 3)
ggsave("10_anno_tsne.pdf", sc_anno1, width = 10, height = 8)
ggsave("10_anno_tsne.png", sc_anno1, width = 10, height = 8)


sc_dot1 <- DotPlot(sc_dat2, features = marker)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
sc_dot1
ggsave("11_dotplot_anno.pdf", sc_dot1, width = 10, height = 8)
ggsave("11_dotplot_anno.png", sc_dot1, width = 10, height = 8)
# sc_type <- DimPlot(luad_data3, reduction = "umap", group.by = "project")

save(sc_dat2, file = "LIHC_sc_anno.rda")

#计算正常样本和疾病样本中不同细胞亚型的比例
# prop.table()
table(sc_dat2$type)
table(Idents(sc_dat2),sc_dat2$type)
cellratio1 <- prop.table(table(Idents(sc_dat2),sc_dat2$type),margin = 2)
cellratio1 <- as.data.frame(cellratio1)
colorcount1 <- length(unique(cellratio1$Var1))
library(ggplot2)
colnames(cellratio1)=c("CellType","Sample","Ratio")
mycolors = c("#FF4500", "#FF7F24", "#556B2F", "#006400", "#20854EFF", "#9AFF9A",   
             "#80B1D3", "#0072B5FF", "#BA55D3", "#E5C494", "#FF0000"
)
#"#FFFFB3","#BEBADA",  "#CD1076", "#ADD8E6","#9400D3","#00BFFF","#E18727FF",   

#绘制算样本中每种细胞的比例的条形图
ratio_plot1 <- ggplot(cellratio1,aes(x=Sample,y=Ratio,fill=CellType))+
  geom_bar(stat="identity",position="fill")+
  #geom_text(aes(label = count,size=8),position=position_fill(vjust =0.5),size = 2.7)+
  scale_fill_manual(values=mycolors)+
  labs(x="Sample",y="Ratio")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(size=0.5),
        text=element_text(size=11,color="black"),
        axis.text=element_text(size=11,color="black"),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.title=element_text(size=16,color="black"))
ratio_plot1
ggsave("12_Disease&control_sample_cellratio_plot.pdf", ratio_plot1, width = 8, height = 6)
ggsave("12_Disease&control_sample_cellratio_plot.png", ratio_plot1, width = 8, height = 6)


#绘制在正常样本和疾病样本中细胞比例差异的箱式图，并计算差异显著性
library(ggplot2)
library(ggpubr)
library(tidyverse)
#增加分组信息
table(Idents(sc_dat2),sc_dat2$orig.ident)
cellratio <- prop.table(table(Idents(sc_dat2),sc_dat2$orig.ident),margin = 2)
cellratio <- as.data.frame(cellratio)
#colorcount <- length(unique(cellratio$Var1))
library(ggplot2)
colnames(cellratio)=c("CellType","Sample","Ratio")
cellratio2 <- cellratio
cellratio2$CellType <- as.character(cellratio2$CellType)
cellratio2$Sample <- as.character(cellratio2$Sample)
range(cellratio2$Ratio)
for(i in 1:dim(cellratio2)[1]){
  if(cellratio2[i,2]=="GSM4909281" | cellratio2[i,2]=="GSM4909282" |
     cellratio2[i,2]=="GSM4909283" | cellratio2[i,2]=="GSM4909284"){
    cellratio2[i,4] <- "BRCA"
  }else{
    cellratio2[i,4] <- "Control"
  }
}

table(cellratio2$V4)
colnames(cellratio2)[4] <- "Group"
# cellratio2$Ratio <- abs(log2(cellratio2$Ratio+1))
range(cellratio2$Ratio)

sample_ratio <- ggplot(cellratio2,aes(x=CellType, y=Ratio, fill=Group)) +
  scale_fill_manual(values = c("#FF7F24","#009ACD")) + #设置颜色
  ylim(0,0.8)+
  geom_boxplot()+
  stat_compare_means(aes(group = Group),method = 'wilcox.test',label = "p.signif",cex=4)+
  xlab("Cell") +
  ylab("Proportion") +
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
  theme(axis.text.x = element_text(face = "bold", color="black",angle = 70,vjust = 1, hjust = 1 ))+#angle = 70表示横轴标题倾斜70度
  theme(axis.text.y = element_text(face = "bold", color="black"))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size = 16),
        legend.title = element_text(size=12))+
  theme(plot.title = element_text(size=16))
sample_ratio
ggsave("13_tumor&control_cell_cellratio_plot.pdf", sample_ratio, width = 8, height = 6)
ggsave("13_tumor&control_cell_cellratio_plot.png", sample_ratio, width = 8, height = 6)

rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./19_scRNA-seq_inferCNV")){
  dir.create("./19_scRNA-seq_inferCNV")
}
setwd("./19_scRNA-seq_inferCNV")

getwd()


library(Seurat)

load("../05_scRNA-seq/LIHC_sc_before_anno.rda")

sc_dat2 <- sa_sc1
sc_dat2 <- RenameIdents(sc_dat2,
                        '0'='T cells', '1'='NKT cells', '2'='T cells',  
                        '3'='T cells', '4'='Macrophage', '5'='NK cells',
                        '6'='Endothelial', '7'='Malignant cells', '8'='T cells',
                        '9'='T cells', '10'='Macrophage', '11'='Plasma cells',
                        '12'='Macrophage', '13'='Fibroblast', '14'='B cells',
                        '15'='Malignant cells', '16'='T cells', '17'='Epithelial', 
                        '18'='DC','19'='Macrophage', '20'='Malignant cells', 
                        '21'='Malignant cells', '22'='Macrophage', '23'='Malignant cells',
                        '24'='Malignant cells', '25'='NK cells', '26'='Malignant cells',
                        '27'='Malignant cells', '28'='Macrophage', '29'='Malignant cells',
                        '30'='Macrophage', '31'='T cells', '32'='Macrophage',
                        '33'='Malignant cells', '34'='Hepatocytes', '35'='NKT cells',
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

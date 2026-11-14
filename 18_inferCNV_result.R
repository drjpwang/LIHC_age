rm(list = ls())
setwd("/data/LIHC/")
if(! dir.exists("./18_inferCNV_result")){
  dir.create("./18_inferCNV_result")
}
setwd("./18_inferCNV_result")

getwd()

library(Seurat)

infer_CNV_obj <- readRDS('../17_inferCNV/inferCNV/run.final.infercnv_obj')
expr <- infer_CNV_obj@expr.data
dim(expr)
expr[1:4,1:4]
data_cnv <- as.data.frame(expr)
dim(expr)
colnames(data_cnv)[1:100]
rownames(data_cnv)[1:100]

#记得要改一下sce中的列名
load("../05_scRNA-seq/LIHC_sc_anno.rda")
UMAP <- sc_dat2
sce <- UMAP[,Idents(UMAP) %in% c("Hepatocytes","Fibroblast")]
colnames(sce) <- gsub("-","_",colnames(sce))
meta <- sce@meta.data


#要根据参数修改哦，我这里的对照组只有Endothelial cells
if(T){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`Fibroblast(non-malignant)`]
  tmp = tmp1
  #tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-2`]
  #tmp= cbind(tmp1,tmp2)
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  cnv_score_mat <- as.matrix(cnv_score_table)
  
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_mat
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  cnv_score_table_pts[1:4,1:4]
  str(as.data.frame(cnv_score_table_pts[1:4,1:4])) 
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
  colnames(cell_scores_CNV) <- "cnv_score" 
}

#可视化
head(cell_scores_CNV)
score=cell_scores_CNV

write.csv(cell_scores_CNV, "cell_scores_CNV.csv")

head(score)
meta$totalCNV = score[match(colnames(sce),
                            rownames(score)),1]
colnames(meta)
p <- ggplot(meta, aes(x= seurat_clusters, 
                      y=totalCNV, 
                      fill=seurat_clusters)) +
  geom_boxplot();print(p)
table(meta$seurat_clusters,meta$cell_type)
ggsave("01_cluster_totalCNV.png",plot = p, width = 8,height = 6,dpi = 600)
ggsave("01_cluster_totalCNV.pdf",plot = p, width = 8,height = 6)

table(meta$seurat_clusters)
str(meta)
meta$seurat_clusters <- as.character(meta$seurat_clusters)

save(meta, file = "Hepatocytes_Fibroblast_scRNA_data.rda")

length(table(meta$seurat_clusters))
cluster <- as.data.frame(table(meta$seurat_clusters))
cluster <- cluster[-1,] #删除参考细胞
cluster$Var1 <- as.character(cluster$Var1)
str(cluster)
seurat_clusters <- cluster$Var1
#计算与参考细胞的拷贝数的差异显著性
for (i in 1:length(seurat_clusters)) {
  meta1 <- meta[meta$seurat_clusters %in% c("13",seurat_clusters[i]),]
  meta1$cellname <- rownames(meta1)
  meta1 <- meta1[,colnames(meta1) %in% c("cellname","totalCNV","seurat_clusters")]
  # meta1 <- meta1[meta1$totalCNV>0,]
  
  sample_ratio <- ggplot(meta1,aes(x=seurat_clusters, y=totalCNV, fill=seurat_clusters)) +
    scale_fill_manual(values = c("#FF7F24","#20B2AA")) + #设置颜色
    # ylim(0,0.8)+
    geom_boxplot()+
    stat_compare_means(aes(group = seurat_clusters),method = 'wilcox.test',label = "p.signif",cex=4)+
    xlab("Clusters") +
    ylab("CNV Score") +
    theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+#去除背景加上黑色外边框
    theme(axis.text.x = element_text(face = "bold", color="black",angle = 70,vjust = 1, hjust = 1 ))+#angle = 70表示横轴标题倾斜70度
    theme(axis.text.y = element_text(face = "bold", color="black"))+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size=14),
          legend.text = element_text(size = 16),
          legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16))
  
  ggsave(paste0("0",i,"_cluster",seurat_clusters[i],"_CNV Score.pdf"),sample_ratio, width = 8,height = 6)
  ggsave(paste0("0",i,"_cluster",seurat_clusters[i],"_CNV Score.png"),sample_ratio, width = 8,height = 6)
  
}

#提取恶性细胞
epi <- UMAP[,Idents(UMAP) %in% c("Hepatocytes")]
#提取表达矩阵
epicount1 <- epi@assays$RNA$count
dim(epicount1)
epicount1 <- as.data.frame(epicount1)
dim(epicount1)
#对所有的恶性细胞进行高低风险组的划分
load("../03_Survival/LASSO_GeneID.rda")
lasso <- data.frame(gene=lasso_geneids$symbol,coefficient=c(0.059247705,0.185107113,0.110993140,
                                                            0.044888230,-0.07235082,0.016215512,
                                                            -0.140727835,-0.035156314))
lasso$coef <- as.numeric(lasso$coefficient)
str(lasso$coef)
rt2 <- epicount1[lasso$gene,]
rt2 <- t(rt2) %>% as.data.frame()
cnt<-1
while (cnt < nrow(rt2)+1) {
  rt2$riskScore[cnt]<-sum(lasso$coef*subset(rt2,select=lasso$gene)[cnt,])
  cnt = cnt + 1
}
str(rt2)
table(is.na(rt2$riskScore))


# 分组 ----------------------------------------------------------------------
# 中位数分组 -------------------------------------------------------------------
rt2 <- rt2[order(rt2$riskScore,decreasing = T),]
dim(rt2)
table(rt2$riskScore)
rt2[(12417+1)/2,]
rt2$risk <- ""
for (i in 1:nrow(rt2)) {
  if(i<(12417+1)/2){
    rt2[i,ncol(rt2)] <- "Low risk"
  }else(
    rt2[i,ncol(rt2)] <- "High risk"
  )
}
# rt2$risk <- ifelse(rt2$riskscore > rt2[10431,3],"High risk","Low risk")
str(rt2)
table(rt2$risk)

write.csv(rt2,"tumor_cell_riskScore.csv")

dim(epi)
risk <- rt2[colnames(epi),]
identical(rownames(risk), colnames(epi))
epi@meta.data$risk <- risk$risk
head(epi@meta.data)


epi1 <- epi
Idents(epi1) <- epi@meta.data$risk
head(epi1@meta.data)
sc_anno <- DimPlot(epi1, reduction = "umap", label = TRUE, pt.size = 5) +
  theme(text=element_text(size= 13 ))#+NoLegend()
ggsave("02_risk_anno_umap.pdf", sc_anno, width = 14, height = 10)
ggsave("02_risk_anno_umap.png", sc_anno, width = 14, height = 10)

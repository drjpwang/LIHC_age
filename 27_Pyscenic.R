rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./27_Pyscenic")){
  dir.create("./27_Pyscenic")
}
setwd("./27_Pyscenic")

getwd()


load("../21_Endothelial/Endothelial_sc_anno.rda")
table(sc_dat2$cell_type)
#使用count
epicount1 <- GetAssayData(sc_dat2, slot = "count", assay = "RNA")
dim(epicount1)
dim(sc_dat2)
df <- data.frame(t(epicount1), check.rows = F, check.names = F)

write.csv(df, "/data/nas1/software/pyscenic/expr_all.csv")


library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)

#### 1.提取 out_SCENIC.loom 信息
loom <- open_loom('/data/nas1/software/pyscenic/aucell.loom')

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

#### 2.加载SeuratData
# UMAP <- get(load('../06_DEC/UMAP.Rdata'))
# UMAP <- get(load('../12_Astrocytes/scRNA.Rdata'))
# UMAP <- readRDS("../11_subtype/keycell_subtype.rds")
# UMAP <- readRDS("../11_pseudotime/FIB_pseudotime.RDS")
# UMAP<- JoinLayers(UMAP)

#选取关注细胞
# key.cells <- subset(UMAP, subset = celltype == "M2_Macrophages")
key.cells <- sc_dat2

sub_regulonAUC <- regulonAUC[,match(colnames(key.cells),colnames(regulonAUC))]
dim(sub_regulonAUC)

#确认是否一致
identical(colnames(sub_regulonAUC), colnames(key.cells))

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression:
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
tf_heatplot <- pheatmap(regulonActivity_byGroup_Scaled,
                        show_rownames = FALSE,  # 显示行名（样本名）
                        show_colnames = TRUE,  # 显示列名（基因名）
                        treeheight_row = 15,   # 行聚类树高度
                        treeheight_col = 15,   # 列聚类树高度
                        color = colorRampPalette(c("blue", "white", "red"))(100),  # 颜色梯度
                        border_color = NA,     # 去除单元格边框
                        fontsize = 20,         # 字体大小
                        main = "TF Activity") # 热图标题

png("01.celltype_TF_plot.png", width = 18, height = 15, units = "in", res = 300)
print(tf_heatplot)
dev.off()

pdf("01.celltype_TF_plot.pdf", width = 18, height = 15)
print(tf_heatplot)
dev.off()

scaled_data=regulonActivity_byGroup_Scaled
head(scaled_data)
library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(scaled_data), function(i){
               dat= data.frame(
                 path  = rownames(scaled_data),
                 cluster =   colnames(scaled_data)[i],
                 sd.1 = scaled_data[,i],
                 sd.2 = apply(scaled_data[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top_gene <- df %>% group_by(cluster) %>% top_n(15, fc)
rowcn = data.frame(path = top_gene$cluster)
scaled_data <- as.data.frame(scaled_data)
n = scaled_data[top_gene$path,]

# 定义颜色标度的范围（假设数据范围是 -2 到 2）
breaks <- seq(-2, 2, length.out = 100)

# 自定义图例标签
tf_heatplot1 <- pheatmap(
  n,
  fontsize = 10,         # 字体大小
  annotation_row = rowcn,
  show_rownames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1),
  breaks = breaks,
  legend_breaks = c(-2, 0, 2)
)

png("02.Artery-1 ECs_TF_plot.png", width = 18, height = 15, units = "in", res = 300)
print(tf_heatplot1)
dev.off()

pdf("02.Artery-1 ECs_TF_plot.pdf", width = 18, height = 15)
print(tf_heatplot1)
dev.off()

endo.markers <- FindAllMarkers(sc_dat2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# endo_Htop100 <- endo.markers %>% group_by(cluster) %>% top_n(n = 150, wt = avg_log2FC) ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少
AE1 <- endo.markers[endo.markers$cluster == "Artery-1 ECs",]
marker <- AE1$gene

AE1_tf <- top_gene[top_gene$cluster=="Artery-1 ECs",]
AE1_tf$gene <- data.frame(sapply(AE1_tf$path, function(x)unlist(strsplit(x,"\\("))[1]),
                          stringsAsFactors=F)[,1]
key_tf <- intersect(marker, AE1_tf$gene)

FeaturePlot(sc_dat2,
            features = key_tf,repel = T, label = F,
            slot = "data", cols = c("lightgrey","red"), ncol = 2, order = T)
dotplot <- DotPlot(sc_dat2, features = key_tf)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
dotplot
ggsave("03.TF_dotplot.pdf", dotplot, width = 10, height = 8)
ggsave("03.TF_dotplot.png", dotplot, width = 10, height = 8)


### 3. rss查看特异TF
selectedResolution <- "celltype"
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
# rssPlot <- plotRSS(rss)
# plotly::ggplotly(rssPlot$plot)
write.csv(rss, "04.rss.csv",  quote = F)
# rss$average_score <- rowMeans(rss, na.rm = TRUE)

rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)


rss <- read.csv("04.rss.csv", row.names = 1)
colnames(rss)
colnames(rss) <- gsub("\\."," ",colnames(rss))
colnames(rss)[2] <- "Artery-2 ECs"
colnames(rss)[4] <- "Artery-1 ECs"
# rss1 <- rss[order(rss$`Artery-1 ECs`,decreasing = T),]
#根据细胞类型绘制热图


# 构建细胞-TF-基因总表
tf_selected <- data.frame()  # 初始化空数据框
tf_gene <- data.frame()


#得到筛选的TF对应的靶基因
key_tf1 <- AE1_tf[AE1_tf$gene %in% key_tf,]
rss1 <- rss[rownames(rss) %in% key_tf1$path, ]
for (i in names(rss1)) {
  rss_order <- rss1[order(-rss1[[i]]), ]
  # rss_order <- rss_order[1:2, ] #TOP5 TF
  tf <- rownames(rss_order)
  for (j in tf) {
    tmp_tf_gene <- data.frame(
      TF = j,
      gene = regulons[j][[j]],
      celltype = i
    )
    tf_gene <- rbind(tf_gene, tmp_tf_gene)
  }
  
  tmp_df <- data.frame(
    TF = rownames(rss_order),
    CellType = i,
    stringsAsFactors = FALSE
  )
  
  tf_selected <- rbind(tf_selected, tmp_df)
}

write.csv(tf_gene, file = "05.celltype_tf_gene.csv", row.names = FALSE, quote = FALSE)

### 4.关键调节子的表达图
# tf_selected <- c("Bhlhe41","Irf2","Nfkb1","Spi1","Yy1")
# tf_selected <- unique(tf_selected$TF)
tf_selected$TF <- gsub("\\(\\+\\)", "", tf_selected$TF)

tf_selected_clean <- unique(tf_selected$TF)


# # 子图
# tf_celltype_map <- tf_selected %>%
#   dplyr::select(TF, CellType) %>%
#   distinct() %>%
#   dplyr::group_by(TF) %>%
#   dplyr::slice(1) %>%
#   tibble::deframe()

# plots <- lapply(tf_selected_clean, function(gene) {
#   celltype <- tf_celltype_map[[gene]]
#   if (is.null(celltype)) celltype <- "Unknown"
#   
#   FeaturePlot(key.cells, features = gene, reduction = "tsne", cols = c("lightgrey", "blue")) +
#     ggtitle(paste0(gene, " (", celltype, ")"))
# })


# 合并图像到一张图
# library(cowplot)
# combined_plot <- plot_grid(plotlist = plots, ncol = 2)

# 保存为 PDF
# pdf("03.tf_selecteds_key.cells_plot.pdf", width = 18, height = 15)
# print(combined_plot)
# dev.off()

# 保存为 PNG
# png("03.tf_selecteds_key.cells_plot.png", width = 18, height = 15, units = "in", res = 300)
# print(combined_plot)
# dev.off()

#网络文件构建
network1 <- tf_gene[,c("TF","celltype")] %>% unique()
network1 <- network1[network1$celltype == "Artery-1 ECs",]
network2 <- tf_gene[,c("TF", "gene")] %>% unique()


write.csv(network1, "06.network_tf_Artery-1 ECs.csv", row.names = F, quote = F)
write.csv(network2, "07.network_tf_gene.csv", row.names = F, quote = F)


#TF和靶基因的通路富集
library(clusterProfiler)
library(org.Hs.eg.db)

#读取候选基因
network3 <- network2
network3$TF1 <- data.frame(sapply(network3$TF, function(x)unlist(strsplit(x,"\\("))[1]),
                           stringsAsFactors=F)[,1]
can_g <- c(network3$gene, network3$TF1)
can_g <- unique(can_g)

#基因ID的转化
pea_gene <- bitr(can_g, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
length(intersect(can_g, pea_gene[,1]))

#GO分析
pea_go <- enrichGO(gene = pea_gene$ENTREZID, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T
)

length(which(pea_go@result[["p.adjust"]]<0.05))
table(pea_go@result[["ONTOLOGY"]])
# BP  CC 
# 1  6

saveRDS(pea_go, "08.candidate_genes_GO.RDS")

write.csv(pea_go, "09.pea_GO.csv")

go_result <- pea_go

#绘制圈图
# 构建GO数据框
go_data = data.frame(Category = go_result$ONTOLOGY,
                     ID = go_result$ID,
                     Term = go_result$Description,
                     Genes = gsub("/", ", ", go_result$geneID),
                     pvalue = go_result$pvalue)
go_data <- go_data[go_data$pvalue<0.05,]  # 7
go_data1 <- go_data[go_data$Category=="BP",]
go_data2 <- go_data[go_data$Category=="CC",]
# go_data3 <- go_data[go_data$Category=="MF",]

go_choose <- rbind(go_data1, go_data2)
# go_choose <- rbind(go_choose, go_data3[1:5,])

#读取差异基因表 
# diffSig <- read.table("../01_GEO_DEG/02_train_diffSig.csv", header = T, row.names = 1, sep = ",", check.names = F)
# id.fc <- diffSig
# geneList <- data.frame(ID = rownames(id.fc), logFC = id.fc$logFC)
# rownames(geneList)=geneList[,1]

library(GOplot)
# circ_GO <- GOplot::circle_dat(go_choose, geneList)
# 
# go_circplot <- GOCircle(data=circ_GO,
#                         nsub=15,
#                         rad1 = 2, rad2 = 3,
#                         zsc.col=c("#96C37D",'white','#5F97D2'),
#                         lfc.col = c('#96C37D','#5F97D2'),
#                         label.size = 4,
#                         label.fontface='bold',
#                         table.legend = T
# )
# ggsave("06_go_circ_plot.pdf", go_circplot, width = 10, height = 8)
# ggsave("06_go_circ_plot.png", go_circplot, width = 10, height = 8)

go_circplot <- clusterProfiler::dotplot(pea_go, split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale="free")
go_circplot

#挑选每个分类中前10显著的条目
BP_top10 <- pea_go %>%
  filter(ONTOLOGY == "BP") %>%  # 筛选 ONTOLOGY 为 BP 的子集
  arrange(pvalue) %>%        # 按照 p.adjust, pvalue 的升序排列
  head(10)                      # 提取前 10行

CC_top10 <- pea_go %>%
  filter(ONTOLOGY == "CC") %>%  # 筛选 ONTOLOGY 为 CC 的子集
  arrange(pvalue) %>%        # 按照 p.adjust 的升序排列
  head(10)                      # 提取前 10行

MF_top10 <- pea_go %>%
  filter(ONTOLOGY == "MF") %>%  # 筛选 ONTOLOGY 为 MF 的子集
  arrange(pvalue) %>%        # 按照 p.adjust 的升序排列
  head(10)                      # 提取前 10行

merge_data <- rbind(BP_top10,CC_top10,MF_top10)#合并三个条目
merge_data$ONTOLOGY <- factor(merge_data$ONTOLOGY, levels = c("BP", "CC","MF"))
merge_data$logPvalue <- (-log(merge_data$pvalue))#数据处理

colnames(merge_data)#检查列名
merge_data$Description <- factor(merge_data$Description, levels = unique(merge_data$Description[order(merge_data$ONTOLOGY)]))

go_barplot <- ggplot(merge_data, aes(x = Description, y = logPvalue, fill = ONTOLOGY)) + #指定数据为merge_data，以Description作为x轴，以logPvalue作为y轴，以ONTOLOGY作为颜色填充。
  geom_bar(stat = "identity", width = 0.8, color = "black") + #添加一个柱状图层，使用identity统计方法（不进行任何变换），每个柱子的宽度为0.8，颜色为黑色
  scale_x_discrete(limits = unique(merge_data$Description[order(merge_data$ONTOLOGY)])) + #设置x轴的标签顺序，根据ONTOLOGY对Description排序，以便按照ONTOLOGY的顺序显示。unique()函数保证每个标签只出现一次。
  coord_flip() + #翻转坐标轴，交换x轴和y轴，使得柱状图变为横向。
  scale_y_continuous(expand = c(0,0))+ #设置y轴的范围为连续型变量（continuous），设置expand参数为c(0,0)，意味着y轴不需要扩展（不需要增加额外的空白）。
  theme(panel.background = element_rect(fill = "white"))+ #设置背景为白色，使用element_rect()函数设置元素的矩形形状
  theme(axis.text = element_text(size = 16, face = "bold", family = "Times"),
        axis.title = element_text(size=18, face = "bold", family="Times",colour = "black"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        legend.title = element_text(size = 18, face = "bold", family = "Times"),
        legend.position = "top")+ #设置各种文本元素的字体、大小和颜色，包括坐标轴文本、坐标轴标题、图例文本和图例标题。将图例放置在顶部。
  labs(x = "GO term", y = "-log P-value")+
  ggsci::scale_fill_nejm(alpha = 0.8) #使用NEJM（New England Journal of Medicine）颜色调色板，设置颜色透明度为0.8，对填充颜色进行比例缩放。
go_barplot
ggsave("06_go_bar_plot.pdf", go_barplot, width = 10, height = 8)
ggsave("06_go_bar_plot.png", go_barplot, width = 10, height = 8)


#KEGG分析
pea_kegg <- enrichKEGG(gene= pea_gene$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
# pea_dotp <- dotplot(pea_kegg, showCategory=10, title="KEGG_PEA_biological")+
#   theme(axis.text.x = element_text(size = 16),
#         axis.text.y = element_text(size = 18),
#         axis.title = element_text(size = 20),
#         plot.title = element_text(size = 20))
length(which(pea_kegg@result[["pvalue"]]<0.05)) # 38个通路显著

saveRDS(pea_kegg, "10.candidate_genes_KEGG.RDS")

# pea_kegg <- readRDS("04_candidate_genes_KEGG.RDS")

#绘制圈图
pea_kegg1 <- setReadable(pea_kegg, 
                         OrgDb = "org.Hs.eg.db", 
                         keyType = "ENTREZID")
which(pea_kegg1@result[["p.adjust"]]<0.05)

kegg_result <- pea_kegg1@result
kegg_result <- as.data.frame(kegg_result)
rownames(kegg_result) <- 1:nrow(kegg_result)
kegg_result <- kegg_result[kegg_result$pvalue < 0.05,] # 22
kegg_result <- kegg_result[order(kegg_result$Count,decreasing = T),]
dim(kegg_result)
write.csv(kegg_result,file = "11.pea_KEGG.csv")


# 气泡图----------------------------------------------------------------------
#根据基因count数排序
dotplot(pea_kegg)
#根据P值排序
df_sorted <- kegg_result[1:10, ]#前10名
colnames(df_sorted)
df_sorted$Description
df_sorted$Description <- gsub(" - Homo sapiens \\(human\\)", "", df_sorted$Description)

kegg_dotplot <- ggplot(df_sorted, aes(x=GeneRatio, y=Description, size=-log10(pvalue), color=pvalue)) + ## 生成散点图，以GeneRatio为x轴，Description为y轴，-log10(p.adjust)为点的大小，p.adjust为点的颜色
  geom_point(alpha=0.7) + # 绘制散点图层，设置透明度为0.7
  scale_color_gradient(low="blue", high="red") +
  labs(title="The most enrichment KEGG pathway", x="GeneRatio", y="KEGG pathway", size="-log10(P-value)", color="P-value") + #设置颜色映射为渐变色，低值为蓝色，高值为红色
  theme_minimal() + #设置图表的标题、x轴和y轴标签，以及点的大小和颜色的图例标签
  theme(axis.text.y = element_text(face = 'bold',size = 12, family = "Times"),
        axis.title = element_text(size = 15,face = 'bold',family = "Times"),
        plot.title = element_text(size = 18, face = "bold", family = "Times")) #设置主题为最小化主题，调整y轴的文本大小和字体粗细，以及标题、轴标签的字体大小
kegg_dotplot
ggsave("12.kegg_dot_plot.pdf", kegg_dotplot, width = 10, height = 8)
ggsave("12.kegg_dot_plot.png", kegg_dotplot, width = 10, height = 8)


network3 <- read.csv("07.network_tf_gene.csv", header = T)
network3$TF1 <- data.frame(sapply(network3$TF, function(x)unlist(strsplit(x,"\\("))[1]),
                           stringsAsFactors=F)[,1]
network3 <- network3[,c(2,3)]
colnames(network3) <- c("Target","TF")
write.csv(network3, "13_TF_Target_network.csv", quote = F, row.names = F)



df=read.table("../00.data/TCGA/TCGA-LIHC.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(4,3,2)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]
##--mRNA--##
gene=unique(network3$TF)
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[gene,]
colnames(data)=substring(colnames(data),1,12)
data=t(data)
sameSample=intersect(row.names(data),row.names(clinical))##365
data=data[sameSample,]
clinical=clinical[sameSample,]
out=cbind(clinical,data)
out=cbind(id=row.names(out),out)
colnames(out)[1] <- "sample"


df_merge <- out

library(survival)
diff_expr_clinical <- df_merge[,-1]
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum


diff_expr_clinical1 <- diff_expr_clinical
colnames(diff_expr_clinical1)
#转录因子
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "ELF1",  #需要计算的数据列名
                                   minprop  = 0.3
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res

diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"ELF1"] > res["ELF1",1],'High','Low')
table(diff_expr_clinical1$risk1)

library(survminer)
# K-M曲线 ----------
kmfit <- survfit(Surv(OS.time, OS) ~ risk1, data = diff_expr_clinical1)
surv_pvalue(kmfit)

pdf('14_KM_ELF1.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High","Low" ),
                     legend.title="Risk",
                     title="ELF1",
                     risk.table = T,
                     risk.table.col = "strata",
                     #linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table +
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "Time (days)") +
  labs(y = "OS survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "bl  ack", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()


png('14_KM_ELF1.png',w=7,h=7,family='Times',units = 'in',res = 600)
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs = c("High","Low" ),
                     legend.title = "Risk",
                     title = "ELF1",
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     # linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette = c("#dc143c","#6495ed"))
s_surv
s_surv$table <- s_surv$table +
  labs(x = "OS Time (Months)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "OS (Months)") +
  labs(y = "OS probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(0,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()


colnames(diff_expr_clinical1)
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "TCF7L1",  #需要计算的数据列名
                                   minprop  = 0.3
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res

diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"TCF7L1"] > res["TCF7L1",1],'High','Low')
table(diff_expr_clinical1$risk1)

library(survminer)
# K-M曲线 ----------
kmfit <- survfit(Surv(OS.time, OS) ~ risk1, data = diff_expr_clinical1)
surv_pvalue(kmfit)

pdf('15_KM_TCF7L1.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High","Low" ),
                     legend.title="Risk",
                     title="TCF7L1",
                     risk.table = T,
                     risk.table.col = "strata",
                     #linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table +
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "Time (days)") +
  labs(y = "OS survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "bl  ack", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()


png('15_KM_TCF7L1.png',w=7,h=7,family='Times',units = 'in',res = 600)
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs = c("High","Low" ),
                     legend.title = "Risk",
                     title = "TCF7L1",
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     # linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette = c("#dc143c","#6495ed"))
s_surv
s_surv$table <- s_surv$table +
  labs(x = "OS Time (Months)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "OS (Months)") +
  labs(y = "OS probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(0,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()



colnames(diff_expr_clinical1)
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "ELK4",  #需要计算的数据列名
                                   minprop  = 0.3
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res

diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"ELK4"] > res["ELK4",1],'High','Low')
table(diff_expr_clinical1$risk1)

library(survminer)
# K-M曲线 ----------
kmfit <- survfit(Surv(OS.time, OS) ~ risk1, data = diff_expr_clinical1)
surv_pvalue(kmfit)

pdf('16_KM_ELK4.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High","Low" ),
                     legend.title="Risk",
                     title="ELK4",
                     risk.table = T,
                     risk.table.col = "strata",
                     #linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table +
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "Time (days)") +
  labs(y = "OS survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "bl  ack", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()


png('16_KM_ELK4.png',w=7,h=7,family='Times',units = 'in',res = 600)
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs = c("High","Low" ),
                     legend.title = "Risk",
                     title = "ELK4",
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     # linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette = c("#dc143c","#6495ed"))
s_surv
s_surv$table <- s_surv$table +
  labs(x = "OS Time (Months)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "OS (Months)") +
  labs(y = "OS probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(0,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()



#转录因子和预后基因的相关性
hubgene <- c("IGFBP3","ELF1","TCF7L1")

df=read.table("../00.data/TCGA/TCGA-LIHC.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(4,3,2)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]

data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[hubgene,]

colnames(data)=substring(colnames(data),1,12)
data=t(data)
sameSample=intersect(row.names(data),row.names(clinical))##365
data=data[sameSample,]

data1 <- data[,c(1,2)]

p1<-ggplot(data1,aes(x=IGFBP3,y=ELF1))+geom_point(aes(x=IGFBP3,y=ELF1), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=data1, method = "spearman",digits = 2)+
  scale_fill_manual(values =c("black"))+
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="IGFBP3",y="ELF1")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1

pdf('17.IGFBP3_ELF1_Correlation.pdf',width=7,height=6,family='Times')
print(p1)
dev.off()

png('17.IGFBP3_ELF1_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p1)
dev.off()# stat.test <- corr.test(res3,method = 'spearman')



data2 <- data[,c(1,3)]

p2<-ggplot(data2,aes(x=IGFBP3,y=TCF7L1))+geom_point(aes(x=IGFBP3,y=TCF7L1), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=data2, method = "spearman",digits = 2)+
  scale_fill_manual(values =c("black"))+
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="IGFBP3",y="TCF7L1")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2

pdf('18.IGFBP3_TCF7L1_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('18.IGFBP3_TCF7L1_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()# stat.test <- corr.test(res3,method = 'spearman')

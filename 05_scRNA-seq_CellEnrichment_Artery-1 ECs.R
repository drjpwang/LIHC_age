rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./05_scRNA-seq_CellEnrichment")){
  dir.create("./05_scRNA-seqCellEnrichment")
}
setwd("./05_scRNA-seqCellEnrichment")

getwd()



if(! dir.exists("./Artery-1 ECs")){
  dir.create("./Artery-1 ECs")
}
setwd("./Artery-1 ECs")

getwd()

#关注细胞的通路富集
load("../../21_Endothelial/Endothelial_sc_anno.rda")

endo.markers <- FindAllMarkers(sc_dat2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(endo.markers, "01_endo.markers.csv", quote = F, row.names = F)

a <- "IGFBP3"

endo_Htop100 <- endo.markers %>% group_by(cluster) %>% top_n(n = 150, wt = avg_log2FC) ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少
AE1 <- endo_Htop100[endo_Htop100$cluster == "Artery-1 ECs",]
setdiff(a, AE1$gene)  # 排名18

marker <- AE1$gene
cell <- c(rep("Artery-1 ECs",nrow(AE1)))

marker_gene <- data.frame(marker=marker, celltype=cell)


data=read.table("../../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)

load("../../03_Survival/train_riskScore.rda")

train_t <- data[, rownames(risk)]


#计算细胞的ssgsea评分
library(GSVA)
gene_list <- split(as.matrix(marker_gene)[,1], marker_gene[,2])
ssgsea_score1 <- gsva(as.matrix(train_t),  gene_list,
                      method = "zscore", kcdf='Gaussian',abs.ranking=TRUE) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

write.csv(ssgsea_score1, "02_Cell_ssgsea_score.csv")

#读取生存信息
df=read.table("../../00.data/TCGA/TCGA-LIHC.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(4,3,2)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]


ssgsea_score1 <- t(ssgsea_score1) %>% as.data.frame()
ssgsea_score2 <- ssgsea_score1
ssgsea_score2$sample <- rownames(ssgsea_score2)

ssgsea_score2$sample <- gsub("\\.", "_", ssgsea_score2$sample)
ssgsea_score2$sample <- gsub("-", "_", ssgsea_score2$sample)

surv <- clinical[rownames(ssgsea_score1),]
surv$sample <- rownames(surv)
surv$sample <- gsub("-", "_", surv$sample)


#整合表达谱数据和生存信息
train_v <- merge(ssgsea_score2, surv, by = "sample")

rownames(train_v) <- train_v$sample
train_v <- train_v[,-1]

# write.csv(surv, file="01_BLCA_survival_infor.csv", quote = F)
# write.csv(train_v, file="02_BLCA_expr_survival.csv", quote = F)

str(train_v)
train_v$OS <- as.numeric(train_v$OS)
train_v$OS.time <- as.numeric(train_v$OS.time)

df_merge <- train_v


# unicox ------------------------------------------------------------------
library(survival)
diff_expr_clinical <- df_merge
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum


diff_expr_clinical1 <- diff_expr_clinical
colnames(diff_expr_clinical1)
#年龄age
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "Artery_1_ECs",  #需要计算的数据列名
                                   minprop  = 0.3
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res # 0.6290213

diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"Artery_1_ECs"] > res["Artery_1_ECs",1],'High','Low')
table(diff_expr_clinical1$risk1)

library(survminer)
library(survival)
# K-M曲线 ----------
range(diff_expr_clinical1$OS.time)
diff_expr_clinical1$OS.time <- diff_expr_clinical1$OS.time/30

kmfit <- survfit(Surv(OS.time, OS) ~ risk1, data = diff_expr_clinical1)
surv_pvalue(kmfit) # GB TW PP mPP FH
pmethod <- surv_pvalue(kmfit, method = "GB")
pmethod


png('03_KM_Artery_1_ECs.png',w=7,h=7,family='Times',units = 'in',res = 600)
s_surv <- ggsurvplot(kmfit,
                     pval = pmethod$pval,   #!!
                     conf.int = F,
                     legend.labs = c("High","Low" ),
                     legend.title = "Risk",
                     title = "OS Curves",
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

pdf('03_KM_Artery_1_ECs.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = pmethod$pval,   #!!
                     conf.int = F,
                     legend.labs = c("High","Low" ),
                     legend.title = "Risk",
                     title = "OS Curves",
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


library(clusterProfiler)
library(org.Hs.eg.db)

#读取候选基因
can_g <- AE1$gene

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
# BP  CC  MF 
# 182  20  10

saveRDS(pea_go, "04_candidate_genes_GO.RDS")

write.csv(pea_go, "05_pea_GO.csv")

go_result <- pea_go

#绘制圈图
# 构建GO数据框
go_data = data.frame(Category = go_result$ONTOLOGY,
                     ID = go_result$ID,
                     Term = go_result$Description,
                     Genes = gsub("/", ", ", go_result$geneID),
                     pvalue = go_result$pvalue)
go_data <- go_data[go_data$pvalue<0.05,]  # 212
go_data1 <- go_data[go_data$Category=="BP",]
go_data2 <- go_data[go_data$Category=="CC",]
go_data3 <- go_data[go_data$Category=="MF",]

go_choose <- rbind(go_data1[1:5,], go_data2[1:5,])
go_choose <- rbind(go_choose, go_data3[1:5,])

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

go_circplot <- clusterProfiler::dotplot(pea_go, split="ONTOLOGY", showCategory =5, colorBy = "pvalue")+
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
  theme(axis.text = element_text(size = 14, face = "bold", family = "Times"),
        axis.title = element_text(size=14, face = "bold", family="Times",colour = "black"),
        legend.text = element_text(size = 14, face = "bold", family = "Times"),
        legend.title = element_text(size = 14, face = "bold", family = "Times"),
        legend.position = "top")+ #设置各种文本元素的字体、大小和颜色，包括坐标轴文本、坐标轴标题、图例文本和图例标题。将图例放置在顶部。
  labs(x = "GO term", y = "-log(P-value)")+
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
length(which(pea_kegg@result[["pvalue"]]<0.05)) # 39个通路显著

saveRDS(pea_kegg, "07_candidate_genes_KEGG.RDS")

# pea_kegg <- readRDS("04_candidate_genes_KEGG.RDS")

#绘制圈图
pea_kegg1 <- setReadable(pea_kegg, 
                         OrgDb = "org.Hs.eg.db", 
                         keyType = "ENTREZID")
which(pea_kegg1@result[["p.adjust"]]<0.05)

kegg_result <- pea_kegg1@result
kegg_result <- as.data.frame(kegg_result)
rownames(kegg_result) <- 1:nrow(kegg_result)
kegg_result <- kegg_result[kegg_result$pvalue < 0.05,] # 39
kegg_result <- kegg_result[order(kegg_result$Count,decreasing = T),]
dim(kegg_result)
write.csv(kegg_result,file = "08_pea_KEGG.csv")


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
  labs(title="KEGG PATHWAY", x="GeneRatio", y="Description", size="-log10(P-value)", color="P-value") + #设置颜色映射为渐变色，低值为蓝色，高值为红色
  theme_minimal() + #设置图表的标题、x轴和y轴标签，以及点的大小和颜色的图例标签
  theme(axis.text.y = element_text(face = 'bold',size = 14, family = "Times"),
        axis.title = element_text(size = 14,face = 'bold',family = "Times"),
        plot.title = element_text(size = 14, face = "bold", family = "Times")) #设置主题为最小化主题，调整y轴的文本大小和字体粗细，以及标题、轴标签的字体大小
kegg_dotplot
ggsave("08_kegg_dot_plot.pdf", kegg_dotplot, width = 10, height = 8)
ggsave("08_kegg_dot_plot.png", kegg_dotplot, width = 10, height = 8)


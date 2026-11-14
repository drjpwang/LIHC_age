rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./14_GSEA_GSVA")){
  dir.create("./14_GSEA_GSVA")
}
setwd("./14_GSEA_GSVA")

getwd()

#GSEA分析
#读取风险得分数据
risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)
risk$sample <- rownames(risk)
risk1 <- risk[,c(ncol(risk)-2,ncol(risk)-1,ncol(risk))]

#读取所有的样本
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-counts.txt",sep="\t",header=T,
                stringsAsFactors = F,quote="",check.names = F,row.names = 1)
range(data)
colnames(data)=substring(colnames(data),1,12)

cnames <- risk1$sample
train_t1 <- data[,cnames]

train_t2 <- as.data.frame(t(train_t1))
train_t2$sample <- rownames(train_t2)

train_t2 <- merge(train_t2, risk1, by = "sample")

#构建分组信息，根据高低风险对样本进行分组
group_t <- data.frame(sample = train_t2$sample, group = train_t2$risk)

rownames(train_t2) <- train_t2$sample
train_t2 <- train_t2[,-c(1,ncol(train_t2)-1,ncol(train_t2))]

table(risk$risk)
table(group_t$group)

#基于高低风险组计算差异基因
train_t2 <- t(train_t2)
dim(train_t2)
library(DESeq2)

#过滤在所有样本都是零表达的基因
train_t2 = train_t2[apply(train_t2,1, function(x) sum(x>1) > 1),]
dim(train_t2)
# 19327   365

dat <- train_t2
dat <- round(dat)
range(dat)
dim(dat)


#设置分组
group_t1 <- group_t
group_t1$group <- factor(group_t1$group, levels = c("High","Low"))
str(group_t1)

table(group_t1$group)

identical(group_t1$sample, colnames(train_t2))

dds<-DESeqDataSetFromMatrix(countData = dat, colData=group_t1, design = ~group)
dds <- dds[rowMeans(counts(dds)) > 1,]
dds <- estimateSizeFactors(dds)
##提取标准化后的数据 
# normalized_counts <- counts(dds,normalized=T)
# 
# write.csv(normalized_counts,file = 'normalized.counts.csv', row.names = T, quote = F)
# save(normalized_counts, file="normalized.counts.rda")

dds<-DESeq(dds)
## 提取差异结果
res =results(dds, contrast = c("group","High","Low"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$pvalue<0.05)
# DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >0.5 )
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff <- 1.5
DEG$regulate <- ifelse(DEG$padj > 0.05, "unchanged",
                       ifelse(DEG$log2FoldChange > logFC_cutoff, "upregulated",
                              ifelse(DEG$log2FoldChange < -logFC_cutoff, "downregulated", "unchanged")))
table(DEG$regulate)
# downregulated     unchanged   upregulated 
#      166         17006           592

#DEG$regulate <- factor(DEG$regulate)

write.csv(DEG, "01_train_DEG.csv")
save(DEG, file="01_train_DEG.rda")

allDiff1 <- DEG

need_DEG <- allDiff1[,c(2,5)] #选择log2FoldChange和pvalue（凑成数据框）

#colnames(need_DEG) <- c('log2FoldChange','pvalue')
need_DEG$SYMBOL <- rownames(need_DEG)

library(GSEABase)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
require(DOSE)
library(Hmisc)
library(future)
library(future.apply)
library(limma)
library(psych)
library(RColorBrewer)
# BiocManager::install("GSVA")
library(GSVA)
# devtools::install_github("junjunlab/GseaVis",force = TRUE)
library(GseaVis)

## 物种设置
organism = 'hsa'    #  人类'hsa' 小鼠'mmu'   
OrgDb = 'org.Hs.eg.db'#人类"org.Hs.eg.db" 小鼠"org.Mm.eg.db"

##### 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）####
#转化id  
df <- bitr(rownames(need_DEG), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = OrgDb)
need_DEG <- merge(need_DEG, df, by='SYMBOL')  #按照SYMBOL合并注释信息
geneList <- need_DEG$log2FoldChange
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)   #从大到小排序
geneList[1:10]

##### gsea富集 ####
KEGG_kk_entrez <- gseKEGG(geneList     = geneList,
                          organism     = organism, #人hsa 鼠mmu
                          pvalueCutoff = 0.25) #实际为padj阈值可调整

KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')#转化id

save(KEGG_kk_entrez, file = "02_survival_GSEA_result.rda")

# GO_kk_entrez <- gseGO(geneList     = geneList,
#                       ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
#                       OrgDb        = OrgDb,#人类org.Hs.eg.db 鼠org.Mm.eg.db
#                       keyType      = "ENTREZID",
#                       pvalueCutoff = 0.25)   #实际为padj阈值可调整
# GO_kk <- DOSE::setReadable(GO_kk_entrez, 
#                            OrgDb=OrgDb,
#                            keyType='ENTREZID')#转化id 



#####  选取富集结果 ####
kk_gse <- KEGG_kk
kk_gse_entrez <- KEGG_kk_entrez

###条件筛选 
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust< 0.25 & abs(kk_gse$NES)>1]
# kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
# kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

#选择展现NES前几个通路 
# down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
# up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$pvalue),decreasing = F),5),]

write.csv(kk_gse_cut,file = "03_survival_gsea_result.csv")

#### 合并 GSEA通路 
diff_gsea$Description
library(enrichplot)
gseap1 <- gseaplot2(kk_gse,
                    diff_gsea$ID,#富集的ID编号
                    title = "Survival_GSEA",#标题
                    color = "red",#GSEA线条颜色
                    base_size = 16,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = F) #显示pvalue等信息
gseap1[[1]] <- gseap1[[1]]+theme(legend.position = "top",legend.direction = "vertical") #将通路显示在图的最上方
gseap1
# ggsave(gseap2, filename = "GSEA_Degs2.pdf",width =15,height =12)

ggexport(gseap1, filename = "04_Survival_GSEA.png",width =800,height =800)
ggexport(gseap1, filename = "04_Survival_GSEA.pdf",width =8,height =8)



#GSVA分析
rm(list = ls())
library(clusterProfiler)
hsgmt <- read.gmt("c2.all.v2023.2.Hs.symbols.gmt") #读gmt文件
gene_list <- split(as.matrix(hsgmt)[,2],
                   hsgmt[,1])

#读取风险得分数据
risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)
risk$sample <- rownames(risk)
risk1 <- risk[,c(ncol(risk)-2,ncol(risk)-1,ncol(risk))]

#读取所有的样本
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)

#提取肿瘤样本
# expr_c <- colnames(train_t)
# for (i in 1:length(expr_c)) {
#   expr_c[i] <- substr(expr_c[i], 1, 15)
#   expr_c[i] <- gsub("\\.", "-", expr_c[i])
# }
# colnames(train_t) <- expr_c

train_t1 <- data

#提取疾病样本
cnames <- intersect(rownames(risk), colnames(train_t1))
train_t1 <- train_t1[,cnames]

group <- risk1[,c("sample","risk")]

table(risk$risk)
table(group$risk)

expr<-train_t1
condition <- group


library(GSVA)
ssgsea.res <-
  gsva(
    as.matrix(expr),
    gene_list,
    method = "ssgsea",
    kcdf = "Gaussian",      
    abs.ranking = T
  )
ssgsea.results<-data.frame(ssgsea.res)
ssgsea.results<-ssgsea.results[order(rownames(ssgsea.results),decreasing = F),]
write.csv(ssgsea.results,'01.gsva.results.csv', quote = FALSE)

expr <- t(expr) %>% data.frame()
expr2 <- cbind(expr,group)

# 分组
library(tidyverse)
library(limma)
group_score <- group$risk %>% as.factor()
design_score <- model.matrix(~0 + group_score)
rownames(design_score) <- rownames(expr)
colnames(design_score) <- levels(group_score)
compare_score <- makeContrasts("High-Low", levels = design_score)

fit <- lmFit(ssgsea.results,design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <-0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0,
         ifelse(allGeneSets$logFC > 0,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$t) > logFCcutoff )
write.csv(allGeneSets,"02.group_c2kegg.csv")
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
write.csv(DEGeneSets,"03.dif_geneset.csv")
##热图
# DEGeneSets <- read.csv('03.dif_geneset.csv', row.names = 1)
# ssgsea.results <- read.csv('01.gsva.results.csv', row.names = 1)


group.gsva<-expr2[order(expr2$risk,decreasing = TRUE),]

rownames(group.gsva) <- str_replace_all(rownames(group.gsva),"-",".")

## 按照分组排一下序
heat<-ssgsea.results[rownames(ssgsea.results)%in%
                       rownames(rbind(head(DEGeneSets[order(DEGeneSets$logFC,decreasing = T),],10),
                                      head(DEGeneSets[order(DEGeneSets$logFC,decreasing = F),],10))),]

heat <- heat[,rownames(group.gsva)]#

group.gsva[1:5,1:5]
heat[1:5,1:5]
# colnames(group.gsva)<-c('group','sample')
# heat <- heat[,group.gsva$id]
# write.csv(rownames(heat),file =paste('0',i,'.',colnames(expr2)[i],'_top20.csv',sep=''),row.names = F)
annotation_col<-as.data.frame(group.gsva$risk)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(heat)
color.key<-c( "#2E8B57", "white", "#CD5C5C")
ann_colors<-list(Group = c(Low="#2E8B57",High="#CD5C5C"))
library(pheatmap)
p<-pheatmap::pheatmap(heat,
                      color = colorRampPalette(color.key)(50),
                      border_color = NA,
                      annotation_col = annotation_col,
                      annotation_colors = ann_colors,
                      labels_row = NULL,
                      clustering_method = 'ward.D2',
                      show_rownames = T,
                      show_colnames = F,
                      fontsize = 6,
                      #fontsize_col = 5,
                      cluster_cols = F,
                      cluster_rows = T,
                      width = 10,
                      #cellwidth=5,
                      # main = paste0(colnames(expr2), " GSVA KEGG")
)
# fn1 = paste0('0',i,'.',colnames(expr2)[i], "_KEGG_heatmap.png")
# fn2 = paste0('0',i,'.',colnames(expr2)[i],"_KEGG_heatmap.pdf")
# png(fn1,w=10,h=6,units = "in",res = 600)
print(p)
# dev.off()
# pdf(fn2,w=10,h=6)
# print(p)
# dev.off()
ggsave("03.hallmaker_gsva.png", p,w=8,h=6)
ggsave("03.hallmaker_gsva.pdf", p,w=8,h=6)



# 发散条形图绘制----------------------------------------------------------------
dat_plot <- data.frame(id = rownames(DEGeneSets),
                       t = DEGeneSets$t)
dat_plot$threshold = factor(DEGeneSets$change)
# dat_plot$threshold = factor(ifelse(dat_plot$p  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
table(dat_plot$threshold)
dat_plot$threshold <- ifelse(dat_plot$threshold == 'UP', 'Up',
                             ifelse(dat_plot$threshold == 'DOWN', 'Down', 'NoSignifi'))
dat_plot <- arrange(dat_plot, t)


dat_plot <-dat_plot[rownames(dat_plot)%in%
                      rownames(rbind(head(dat_plot[order(dat_plot$t,decreasing = T),],10),
                                     head(dat_plot[order(dat_plot$t,decreasing = F),],10))),]


dat_plot <- dat_plot[which(abs(dat_plot$t) > 2), ] # 筛选t阈值大于2的基因集
length(dat_plot$id)
dat_plot$id <- gsub('KEGG_', '', dat_plot$id)
dat_plot$id <- factor(dat_plot$id, levels=dat_plot$id)

library(ggplot2)
library(dplyr)

# 绘制柱状图
library(ggprism)

p <- ggplot(data = dat_plot, aes(x = id, y = t, fill = threshold)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c('Up' = '#FF0088', 'NoSignifi' = '#cccccc', 'Down' = '#7bcd7b')) +
  scale_y_continuous(breaks = seq(-8, 15, 2)) +
  xlab('') +
  ylab('t value of GSVA score') + # 注意坐标轴旋转了
  guides(fill = guide_legend(title = "Threshold")) + # 显示图例并添加标题
  theme_prism(border = TRUE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# 计算标签数量
Sepsis1 <- dat_plot %>% filter(t < -1) %>% nrow()
Sepsis0 <- dat_plot %>% filter(t < 0) %>% nrow()
Sepsis_ARDS0 <- dat_plot %>% filter(t < 1) %>% nrow()
Sepsis_ARDS1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:Sepsis1, ], aes(x = id, y = 0.1, label = id),
                   hjust = 0, color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(Sepsis_ARDS0 + 1):Sepsis_ARDS1, ], aes(x = id, y = -0.1, label = id),
            hjust = 1, color = 'black') # 大于1的为黑色标签

# 保存为PDF
pdf(file = '04_GSVA_bar.pdf', width = 21, height = 10)
print(p)
dev.off()

# 保存为PNG
png(file = '04_GSVA_bar.png', width = 21, height = 10, units = 'in', res = 600)
print(p)
dev.off()

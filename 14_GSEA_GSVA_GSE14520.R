rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./14_GSEA_GSVA_GSE14520")){
  dir.create("./14_GSEA_GSVA_GSE14520")
}
setwd("./14_GSEA_GSVA_GSE14520")

getwd()


#GSVA分析
library(clusterProfiler)
hsgmt <- read.gmt("../14_GSEA_GSVA/c2.all.v2023.2.Hs.symbols.gmt") #读gmt文件
gene_list <- split(as.matrix(hsgmt)[,2],
                   hsgmt[,1])

#读取风险得分数据
risk <- read.csv("../04_GSE14520_SurvivalVerify/verify_riskScore.csv", header = T, row.names = 1)
risk$sample <- rownames(risk)
risk1 <- risk[,c("sample","risk")]

#读取所有的样本
data=read.table("../00.data/GEO_HCC/GSE14520/GSE14520_GPL3921.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)

data1 <- data[,colnames(data) %in% risk1$sample]

#提取肿瘤样本
# expr_c <- colnames(train_t)
# for (i in 1:length(expr_c)) {
#   expr_c[i] <- substr(expr_c[i], 1, 15)
#   expr_c[i] <- gsub("\\.", "-", expr_c[i])
# }
# colnames(train_t) <- expr_c

train_t1 <- data1

#提取疾病样本
identical(colnames(train_t1), risk1$sample)

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

rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./09_ImmuneInfiltration")){
  dir.create("./09_ImmuneInfiltration")
}
setwd("./09_ImmuneInfiltration")

getwd()


#读取所有的样本
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)


#读取风险得分数据
risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)

train_expr <- data[,rownames(risk)]


write.csv(train_expr, file = "TCGA-LIHC-mRNA-tpm.csv")

train_t <- train_expr

source("/data/nas1/software/cibersort/Cibersort_source.R")

# sig_matrix <- read.table("/data/nas1/zhoujie/software/cibersort/LM22.txt", sep = "\t")

#运行CIBERSORT
sig_res <- CIBERSORT("/data/nas1/software/cibersort/LM22.txt", "TCGA-LIHC-mRNA-tpm.csv", 
                     perm = 1000, QN = F)
#QN如果是芯片设置为T，如果是测序就设置为F

write.csv(sig_res, "01_cibersort_res.csv", 
          row.names = T, quote = F)


sig_res[1:3,22:23]
#将CIBERSORT数据整理
sig_res1 <- sig_res
sig_res1 <- sig_res1[,1:22]
sig_res1 <- t(sig_res1)
sig_cell <- as.data.frame(row.names(sig_res1))
colnames(sig_cell) <- "Cell"
sig_res1 <- as.data.frame(sig_res1)
sig_res1 <- cbind(sig_cell, sig_res1)
rownames(sig_res1) <- 1:dim(sig_res1)[1]

#pdf('02.cibersort.box.pdf',w=10,h=8)
# for (i in 1:dim(group_t)[1]) {
#   group_t[i,1] <- substr(group_t[i,1], 1, 15)
#   group_t[i,1] <- gsub("\\.", "-", group_t[i,1])
# }

#根据高低风险组进行分组
risk1 <- data.frame(risk=risk[,ncol(risk)])
#library(dplyr)
library(tidyverse)
rownames(risk1) <- rownames(risk)
risk1$sample <- rownames(risk1)
sig_res1 <- sig_res1 %>%
  gather(sample, fraction, -Cell) %>%
  merge(risk1,by='sample')

write.csv(sig_res1, "02_cibersort_result_cell&group.csv")

#PRAD患者高低风险组样本之间的免疫浸润差异柱状堆积图
library(RColorBrewer)
# mypalette <- colorRampPalette(brewer.pal(7,"Accent"))
mypalette <- c(
  "#E41A1C",  # 红色
  "#377EB8",  # 蓝色
  "#4DAF4A",  # 绿色
  "#984EA3",  # 紫色
  "#FF7F00",  # 橙色
  "#FFFF33",  # 黄色
  "#A65628",  # 棕色
  "#F781BF",  # 粉红色
  "#999999",  # 灰色
  "#66C2A5",  # 浅绿色
  "#FC8D62",  # 浅橙色
  "#8DA0CB",  # 浅蓝色
  "#E78AC3",  # 浅粉红色
  "#A6D854",  # 浅黄色绿色
  "#FFD92F",  # 浅黄色
  "#E5C494",  # 浅棕色
  "#B3B3B3",  # 浅灰色
  "#1F78B4",  # 深蓝色
  "#33A02C",  # 深绿色
  "#FB9A99",  # 浅红色
  "#CAB2D6",  # 浅紫色
  "#6A3D9A"   # 深紫色
)                                            # 高对比度

# 可视化调色板
barplot(rep(1, length(mypalette)), col = mypalette, space = 0, border = NA)
ciber_histogram <- ggplot(sig_res1, aes(x=sample, y=fraction, fill=Cell)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette)+
  facet_grid(~risk,scales= "free",space= "free")
#dev.off()
ciber_histogram
ggsave("03_cibersort_Histogram.pdf", ciber_histogram, width = 10, height = 8)
ggsave("03_cibersort_Histogram.png", ciber_histogram, width = 10, height = 8)

#PRAD患者样本和对照样本之间的免疫浸润差异箱式图
sig_res2 <- sig_res
sig_res2 <- sig_res2[,1:22]
sig_res2 <- t(sig_res2)
sig_cell <- as.data.frame(row.names(sig_res2))
sig_res2 <- sig_res2[which(rowSums(sig_res2)>0),]

sig_cell <- as.data.frame(row.names(sig_res2))
colnames(sig_cell) <- "Cell"

sig_res2 <- cbind(sig_cell, sig_res2)
rownames(sig_res2) <- 1:dim(sig_res2)[1]
write.csv(sig_res2, "04_cibersort_result_cell.csv")

# keep <- rowSums(sig_res2>0)>=floor(0.3*ncol(sig_res2))  # 过滤掉在所有样本中占比低于30%的细胞类型
# sig_res2 <- sig_res2[keep,]

#剩下14种免疫细胞
dat.cibersort <- sig_res2 %>% 
  #tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
colnames(dat.cibersort) <- dat.cibersort[1,]
dat.cibersort <- dat.cibersort[-1,]
colnames(dat.cibersort)[1] <- "sample"

dat.cibersort <- merge(risk1, dat.cibersort, by = "sample")

# dat.cibersort <- subset(dat.cibersort, select=-c(2))
dat.cibersort1 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "risk"))

colnames(dat.cibersort1)
library(rstatix)

dat.cibersort1$Score <- as.numeric(dat.cibersort1$Score)
stat_cibersort <- dat.cibersort1 %>% 
  group_by(ImmuneCell) %>%
  rstatix::wilcox_test(Score ~ risk) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.csv(stat_cibersort,file = '05_cibersort_ImmuneCell.csv',row.names = T)

DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort,file = '06_cibersort_DE.csv',row.names = T)

colnames(dat.cibersort1)
#violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
violin.cibersort<-dat.cibersort1
library(ggpubr)
ciber_boxplot <- ggplot(dat.cibersort1, aes(x=ImmuneCell,
                                            y=Score,
                                            fill=risk)) +
  # geom_violin(trim=F,color="black") + "#FF7F24","#5F9EA0"
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values= c("#dc143c","#6495ed"), name = "Risk")+
  labs(title="Immune Infiltration", x="", y = "Fraction",size=20) +
  stat_compare_means(data = dat.cibersort1,
                     mapping = aes(group = risk),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=60,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+facet_wrap(~ImmuneCell,scales = "free",nrow = 3) 
ciber_boxplot
ggsave("07_cibersort_boxplot.pdf", ciber_boxplot, width = 14, height = 8)
ggsave("07_cibersort_boxplot.png", ciber_boxplot, width = 14, height = 8)

#LUSC中预后基因与差异免疫细胞相关性分析
tiics_result <- read.csv('04_cibersort_result_cell.csv',check.names = F, row.names = 1, header = T)
# %>% lc.tableToNum
#tiics_result <- tiics_result[, group$sample] %>% as.matrix()
rownames(tiics_result) <- tiics_result$Cell
tiics_result <- tiics_result[,-1]
tiics_result <- t(tiics_result) %>% as.data.frame()
# 
# #tiics_result <- tiics_result[rownames(expr),]
# 
diff <- read.csv('06_cibersort_DE.csv', row.names = 1, header = T)
tiics_result <- tiics_result[,diff$ImmuneCell]
tiics_result1 <- tiics_result
env.cor <- round(cor((tiics_result1),method="spearman"), 3)
# env.p <-round(cor_pmat((gene_dat),method = "spearman"),3)
cor_p <- WGCNA::corPvalueStudent(env.cor,nrow(tiics_result1))

library(corrplot)

#PRAD中差异免疫细胞与差异免疫细胞相关性分析
pdf("07_cibersort_ImmuneCell_corplot.pdf", width = 12, height = 8)
corrplot(corr =env.cor,p.mat = cor_p,type="upper",
         col = colorRampPalette(c("blue", "white", "red"))(50),
         tl.pos="lt",tl.col="black",
         insig = "label_sig", sig.level = c(.001,.01, .05),
         pch.cex=1,pch.col = "black",order = "AOE")
corrplot(corr = env.cor,type="lower",add=TRUE,method="number",
         col = colorRampPalette(c("blue", "white", "red"))(50),
         tl.pos="n",tl.col="black",tl.cex=1.2,
         diag=FALSE, cl.pos="n",pch.col = "black",
         number.cex = 1.2,order = "AOE")
dev.off()

png("07_cibersort_ImmuneCell_corplot.png", width = 1200, height = 800)
corrplot(corr =env.cor,p.mat = cor_p,type="upper",
         col = colorRampPalette(c("blue", "white", "red"))(50),
         tl.pos="lt",tl.col="black",
         insig = "label_sig", sig.level = c(.001,.01, .05),
         pch.cex=1,pch.col = "black",order = "AOE")
corrplot(corr = env.cor,type="lower",add=TRUE,method="number",
         col = colorRampPalette(c("blue", "white", "red"))(50),
         tl.pos="n",tl.col="black",tl.cex=1.2,
         diag=FALSE, cl.pos="n",pch.col = "black",
         number.cex = 1.2,order = "AOE")
dev.off()

#BRCA中预后基因与差异免疫细胞相关性分析
library(corrplot)
#读取表达数据
#load("01_train_expression.rda")

#读取预后基因
lasso_geneids <- readRDS("../03_Survival/08_lasso_geneids.rds")
survival_gene <- lasso_geneids

expr <- t(train_t) %>% as.data.frame()
range(expr)
#expr <- log2(expr+1)
expr <- expr[,survival_gene]

tiics_result <- read.csv('04_cibersort_result_cell.csv',check.names = F, row.names = 1, header = T)

rownames(tiics_result) <- tiics_result$Cell
tiics_result <- tiics_result[,-1]
tiics_result <- t(tiics_result) %>% as.data.frame()
#tiics_result <- tiics_result[rownames(expr),]

diff <- read.csv('06_cibersort_DE.csv', row.names = 1, header = T)
tiics_result <- tiics_result[,diff$ImmuneCell]

#比较两个对象是否相等
identical(rownames(expr), rownames(tiics_result))
expr <- expr[rownames(tiics_result),]

cor_r <- cor(expr,tiics_result,method = "spearman")
#cor_p <- WGCNA::corPvalueStudent(cor_r,length(rownames(dat_exp_diff)))
library(psych)
d <- corr.test(expr,tiics_result,use="complete",method = 'spearman')
cor_p1 <- d$p
cor_r1 <- cor_r %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>%
  tidyr::gather(., cell,Correlation,-gene)
cor_p2 <- cor_p1 %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") %>% 
  tidyr::gather(., cell, Pvalue, -gene)

cor_dat <- cbind(cor_r1, cor_p2)[,c("gene","cell","Correlation","Pvalue")]
#cor_dat <- cor_dat[cor_dat$Cell %in% stat_res$Cell,]
write.csv(cor_dat,"08_correlation_cor.csv")


cor_dat <- cor_dat %>%
  mutate(text = case_when(
    Pvalue <= 0.001 ~ "\n***",
    between(Pvalue, 0.001, 0.01) ~ "\n**",
    between(Pvalue, 0.01, 0.05) ~ "\n*",
    T ~ ""))
# data <- data %>%
#   mutate(text = case_when(
#     Pvalue > 0 ~ paste(round(Pvalue, 2), "+"), # round()
#     Pvalue < 0 ~ paste(round(Pvalue, 2), "-")))
ciber_corplot3 <- 
  ggplot(cor_dat, aes(gene, cell)) + 
  geom_tile(aes(fill = Correlation), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") +
  geom_text(aes(label = text),col ="black",size = 5) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold")) +
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) + 
  scale_x_discrete(position = "top")
ciber_corplot3
ggsave("09_cibersort_Biomarker_ImmuneCell_corplot.pdf", ciber_corplot3, width = 14, height = 8)
ggsave("09_cibersort_Biomarker_ImmuneCell_corplot.png", ciber_corplot3, width = 14, height = 8)


r_result <- corr.test(expr,tiics_result,method = 'spearman')

r <-
  r_result$r %>%
  melt() %>%
  magrittr::set_names(c('tax', 'Index', 'r'))

p <-
  r_result$p %>%
  melt() %>%
  magrittr::set_names(c('tax', 'Index', 'P_value')) %>%
  mutate(P_value_sig = case_when(P_value > 0.05 ~ "ns",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))


data <- cbind(r,p) %>% dplyr::select(-(4:5)) #删除重复列
head(data)
# write.csv(data,"05_cor_cell_gene.csv")
# # 散点图
# for (i in 1:4) {
#   # i <- 
#   gene <- data[i,1]
#   CELL <- data[i,2]
#   x <- de.dat[,CELL]
#   y <- hub.exp[,gene]
#   df <- cbind(x,y)%>%as.data.frame()
#   p1 <- ggplot(df,aes(x,y))+
#     xlab(gene)+ylab(CELL)+
#     geom_point()+geom_smooth(method = 'lm',formula = y~x)+
#     theme_bw()+
#     stat_cor(method = 'spearman',aes(x=x,y=y))
#   p1
#   ggsave(p1,filename=paste0(i+2,".",gene,"_",CELL,'.pdf'),w=5,h=5)
#   ggsave(p1,filename=paste0(i+2,".",gene,"_",CELL,'.png'),w=5,h=5)
# }

textMatrix =  paste(round(data$r,3), "\n(",
                    data$P_value_sig, ")", sep = "")
# 绘图主题
# theme <- theme_bw(base_size = 15) +
#   theme(panel.grid.major = element_blank(),
#         text = element_text(face='bold'),
#         legend.key.size = unit(8, "pt"),
#         panel.grid.minor = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.x=element_text(angle=45,vjust=1, hjust=1))

immune_corplot <- ggplot(data,aes(x=Index, y=tax))+
  geom_tile(aes(fill=r), color = 'white',alpha = 0.6) +
  geom_text(label = textMatrix, color = 'black', size =4) +
  ggsci::scale_fill_gsea() +
  theme(legend.position="top") +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(),
        text = element_text(face='bold'),
        legend.key.size = unit(8, "pt"),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(hjust = 1, size = 12, face = "bold", angle=60),
        axis.text.y = element_text(size = 12, face = "bold")) + 
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) + 
  # scale_x_discrete(position = "top") + 
  xlab('') +
  ylab('')
# ggsave('heatmap_corr.png',height = 2, width = 2)
immune_corplot


####优化图片
fill_scale <- ggsci::scale_fill_gsea(
  name = "Spearman Correlation",
  breaks = c(-0.2, 0, 0.2, 0.3), # 匹配原图图例的数值标记
  labels = c("-0.2", "0", "0.2", "0.3")
)


# 3. 绘图（完善图例+美化细节）
immune_corplot <- ggplot(data, aes(x = Index, y = tax)) +
  # 热图色块（调整边框和透明度）
  geom_tile(aes(fill = r), color = "white", alpha = 0.8, linewidth = 0.8) +
  # 文本标签（优化大小和加粗）
  geom_text(
    label = textMatrix, 
    color = "black", 
    size = text_size, 
    fontface = "bold",
    vjust = 0.5, hjust = 0.5
  ) +
  # 配色+图例美化（明确数值范围）
  fill_scale +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 12, barheight = 1.2, # 调整图例宽度高度
      label.position = "bottom",
      ticks = TRUE # 显示刻度线，对应数值标记
    )
  ) +
  # 主题优化（更简洁专业）
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(), # 移除所有网格线
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    # 坐标轴文本优化（避免重叠）
    axis.text.x = element_text(
      angle = 60, hjust = 1, vjust = 1,
      size = 12, face = "bold", color = "black"
    ),
    axis.text.y = element_text(
      size = 12, face = "bold", color = "black"
    ),
    # 图例美化
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "top",
    # 调整边距（避免元素被截断）
    plot.margin = ggplot2::margin(10, 10, 10, 10, unit = "pt")
  ) +
  # 显著性说明（放在图例标题下方）
  labs(
    fill = paste0(
      "Spearman Correlation\n",
      "* p < 0.05  |  ** p < 0.01  |  *** p < 0.001"
    )
  )
immune_corplot
ggsave(filename = '10_Biomarker_ImmuneCell_corplot.pdf',immune_corplot,w=10,h=8)
ggsave(filename = '10_Biomarker_ImmuneCell_corplot.png',immune_corplot,w=10,h=8)

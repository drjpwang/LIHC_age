rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./11_ImmuneInfunction")){
  dir.create("./11_ImmuneInfunction")
}
setwd("./11_ImmuneInfunction")

getwd()

#读取免疫细胞相关基因
immune_c <- read.table("ImmuneCell_geneset.csv", header = T, sep = ",")
colnames(immune_c)[1] <- "GeneName"

#读取免疫功能相关基因
immune_f <- read.table("ImmuneFunctionl_geneset.csv", header = T, sep = ",")
colnames(immune_f) <- c("GeneName", "CellFunction")

#ssgse分析
#读取风险得分数据
risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)
risk$sample <- rownames(risk)
risk1 <- risk[,c("sample","risk")]

#输入疾病样本表达矩阵fpkm
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)

cnames <- risk1$sample
train_t <- data[,cnames]

gene_list1 <- split(as.matrix(immune_c)[,1], immune_c[,2])
gene_list2 <- split(as.matrix(immune_f)[,1], immune_f[,2])

#计算免疫细胞的ssgsea评分
library(GSVA)
ssgsea_score1 <- gsva(as.matrix(train_t),  gene_list1,
                      method = "ssgsea", kcdf='Gaussian',abs.ranking=TRUE) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

write.csv(ssgsea_score1, "01_ImmuneCell_ssgsea_score.csv")

#计算免疫功能的ssgsea评分
ssgsea_score2 <-  gsva(as.matrix(train_t),  gene_list2,
                       method = "ssgsea", kcdf='Gaussian',abs.ranking=TRUE) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

write.csv(ssgsea_score2, "02_ImmuneFunction_ssgsea_score.csv")


res.cibersort1 <- data.frame(ssgsea_score1)
colnames_sum1 <- colnames(ssgsea_score1)
colnames(res.cibersort1) <- colnames_sum1
res.cibersort1$cell_type <- rownames(res.cibersort1)
rownames(res.cibersort1) <- NULL

library(tidyr)
## 高低风险组的差异免疫细胞-------
dat.cibersort1 <- res.cibersort1 %>% 
  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
# group <- data.frame(sample = riskScore$id,group = riskScore$risk)
dat.cibersort1 <- merge(risk1, dat.cibersort1, by = "sample")
###删除所有值均为0的列 purrr::discard(df, ~all(.x == 0)) purrr::keep(df, ~any(.x != 0))
dat.cibersort1 <- dat.cibersort1[,apply(dat.cibersort1,2,function(dat.cibersort1) !all(dat.cibersort1==0))]
dat.cibersort2 <- tidyr::gather(dat.cibersort1, ImmuneCell, Score, -c("sample", "risk"))

colnames(dat.cibersort2)
library(rstatix)
stat_cibersort <- dat.cibersort2 %>% 
  dplyr::group_by(ImmuneCell) %>%
  rstatix::wilcox_test(Score ~ risk) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr   ###"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  add_significance("p")
write.csv(stat_cibersort,file = '03_stat_cell.csv')


DE.cibersort1<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort1,file = '04_DE_cell.csv')

colnames(dat.cibersort2)
violin.cibersort1<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]

violin.cibersort1$risk <- factor(violin.cibersort1$risk,levels =c("High","Low") )

library(ggplot2)
library(ggrepel)
library(ggpubr)
violin_plot1 <- ggplot(violin.cibersort1, aes(x=ImmuneCell, y=Score,fill=risk)) +
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  # geom_violin(aes(fill=group),trim=FALSE,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  # geom_point(aes(color = group), position=position_dodge(0.9), size = 1) +
  scale_fill_manual(values= c("#dc143c","#6495ed"), name = "Risk")+
  labs(title="Immune Cell", x="", y = "",size=20) +
  stat_compare_means(data = violin.cibersort1,
                     mapping = aes(group = risk),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  labs(y = "Immune score")+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot1 # 10种差异免疫细胞

ggsave(filename = '05_cell_violin.pdf',violin_plot1,w=8,h=6)
ggsave(filename = '05_cell_violin.png',violin_plot1,w=8,h=6)
#dev.off()


res.cibersort2 <- data.frame(ssgsea_score2)
colnames_sum2 <- colnames(ssgsea_score2)
colnames(res.cibersort2) <- colnames_sum2
res.cibersort2$function_type <- rownames(res.cibersort2)
rownames(res.cibersort2) <- NULL

## 高低风险组的差异免疫功能-------
dat.cibersort3 <- res.cibersort2 %>% 
  tibble::column_to_rownames(var = "function_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
# group <- data.frame(sample = riskScore$id,group = riskScore$risk)
dat.cibersort3 <- merge(risk1, dat.cibersort3, by = "sample")
###删除所有值均为0的列 purrr::discard(df, ~all(.x == 0)) purrr::keep(df, ~any(.x != 0))
dat.cibersort3 <- dat.cibersort3[,apply(dat.cibersort3,2,function(dat.cibersort3) !all(dat.cibersort3==0))]
dat.cibersort4 <- tidyr::gather(dat.cibersort3, ImmuneFunction, Score, -c("sample", "risk"))

colnames(dat.cibersort4)
library(rstatix)
stat_cibersort <- dat.cibersort4 %>% 
  dplyr::group_by(ImmuneFunction) %>%
  rstatix::wilcox_test(Score ~ risk) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr   ###"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  add_significance("p")
write.csv(stat_cibersort,file = '06_stat_function.csv')


DE.cibersort2<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort2,file = '07_DE_function.csv')

colnames(dat.cibersort4)
violin.cibersort2<-dat.cibersort4[dat.cibersort4$ImmuneFunction%in%stat_cibersort$ImmuneFunction[which(stat_cibersort$p<0.05)],]

violin.cibersort2$risk <- factor(violin.cibersort2$risk,levels =c("High","Low") )

violin_plot2 <- ggplot(violin.cibersort2, aes(x=ImmuneFunction, y=Score,fill=risk)) +
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  # geom_violin(aes(fill=group),trim=FALSE,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  # geom_point(aes(color = group), position=position_dodge(0.9), size = 1) +
  scale_fill_manual(values= c("#dc143c","#6495ed"), name = "Risk")+
  labs(title="Immune Function", x="", y = "",size=20) +
  stat_compare_means(data = violin.cibersort2,
                     mapping = aes(group = risk),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  labs(y = "Immune score")+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot2 # 6种差异免疫功能

ggsave(filename = '08_cell_violin.pdf',violin_plot2,w=8,h=6)
ggsave(filename = '08_cell_violin.png',violin_plot2,w=8,h=6)
#dev.off()

#计算风险评分与差异免疫细胞和差异免疫功能的相关性
riskScore <- data.frame(riskScore = risk1$riskScore)
rownames(riskScore) <- risk1$sample

# ---1--- 差异免疫细胞
ssgsea_score1 <- t(ssgsea_score1) %>% as.data.frame()

library(psych)
corr<- psych::corr.test(riskScore,ssgsea_score1,method =  'spearman',adjust = 'none')
#m=corr$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(de.dat))+1):ncol(nc)]
cor_r<-as.data.frame(corr$r)
#p=corr$p[1:nrow(cor.dat),(ncol(nc)-length(rownames(de.dat))+1):ncol(nc)]
cor_p<-as.data.frame(corr$p)

cor_r2 <- cor_r %>% tibble::rownames_to_column(var = "riskScore") %>%
  tidyr::gather(., cell_type,Correlation,-riskScore)#转换数据长短 cor
cor_p2 <- cor_p  %>% tibble::rownames_to_column(var = "riskScore") %>%
  tidyr::gather(., cell_type, Pvalue, -riskScore)#转换数据长短 p
cor_dat <- cbind(cor_r2, cor_p2)[,c("riskScore","cell_type","Correlation","Pvalue")]

write.csv(cor_dat,file ='09_Correlation.ImmuneCell.csv',row.names = F,quote = F)

#tmp<- cor_dat[cor_dat$Pvalue <0.05 & abs(cor_dat$Correlation)>0.3,]
# tmp<- cor_dat

#棒棒糖图--------------

correlation1 <- cor_dat
# if (nrow(correlation)==0){next}
correlation1$color<- ifelse(abs(correlation1$Correlation) >0.3 & correlation1$Pvalue<0.05,'red','black')
correlation1$label <- ifelse(abs(correlation1$Correlation) >0.3 & correlation1$Pvalue<0.05,sprintf("%0.2e",correlation1$Pvalue),"")
correlation1<- correlation1[order(-correlation1$Correlation),]

p0<-ggdotchart(correlation1, x = "cell_type", y = "Correlation",
               dot.size = 4,
               # filled.contour(p.value,col=Lab.palette(20)),
               color ='Pvalue',
               label = 'label',
               font.label = list(color = correlation1$color,size = 8),
               #ylab="mgp = c(3, 2, 0)",
               # font.label = list(size=10,face="bold",color='black',position="dodge",
               #                   vjust=0.5),
               #  y.label=0.7,
               # palette = palette,
               # 按照cyl填充颜色
               # palette = palette, # 修改颜色
               sorting = "descending",
               
               add = "segments",                             # 添加棒子
               rotate = TRUE,
               ggtheme = theme_pubr(),                        # 改变主题
               ylab="Correlation",
               xlab='',
               # dot.size = 6 ,
               title="RiskScore_ImmuneCell_Correlation"
)+scale_colour_gradient( high = "yellow",low = "#246BAE")
p10<-p0+
  theme(legend.position = "right",
        panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
  theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
        axis.text.x =element_text(size=12,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=12,family = "Times", face = "bold"),
        plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
        legend.text = element_text(size = 16, family = "Times"),
        legend.title = element_text(size = 18, family = "Times",face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p10
ggsave("10_RiskScore_ImmuneCell_Correlation.pdf",height = 6, width = 10, p10)
ggsave("10_RiskScore_ImmuneCell_Correlation.png",height = 6, width = 10, p10)


# ---2--- 差异免疫功能
ssgsea_score2 <- t(ssgsea_score2) %>% as.data.frame()

library(psych)
corr<- psych::corr.test(riskScore,ssgsea_score2,method =  'spearman',adjust = 'none')
#m=corr$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(de.dat))+1):ncol(nc)]
cor_r<-as.data.frame(corr$r)
#p=corr$p[1:nrow(cor.dat),(ncol(nc)-length(rownames(de.dat))+1):ncol(nc)]
cor_p<-as.data.frame(corr$p)

cor_r2 <- cor_r %>% tibble::rownames_to_column(var = "riskScore") %>%
  tidyr::gather(., cell_function,Correlation,-riskScore)#转换数据长短 cor
cor_p2 <- cor_p  %>% tibble::rownames_to_column(var = "riskScore") %>%
  tidyr::gather(., cell_function, Pvalue, -riskScore)#转换数据长短 p
cor_dat <- cbind(cor_r2, cor_p2)[,c("riskScore","cell_function","Correlation","Pvalue")]

write.csv(cor_dat,file ='11_Correlation.ImmuneFunction.csv',row.names = F,quote = F)

#tmp<- cor_dat[cor_dat$Pvalue <0.05 & abs(cor_dat$Correlation)>0.3,]
# tmp<- cor_dat

#棒棒糖图--------------

correlation2 <- cor_dat
# if (nrow(correlation)==0){next}
correlation2$color<- ifelse(abs(correlation2$Correlation) >0.3 & correlation2$Pvalue<0.05,'red','black')
correlation2$label <- ifelse(abs(correlation2$Correlation) >0.3 & correlation2$Pvalue<0.05,sprintf("%0.2e",correlation1$Pvalue),"")
correlation2<- correlation2[order(-correlation2$Correlation),]

p0<-ggdotchart(correlation2, x = "cell_function", y = "Correlation",
               dot.size = 4,
               # filled.contour(p.value,col=Lab.palette(20)),
               color ='Pvalue',
               label = 'label',
               font.label = list(color = correlation2$color,size = 8),
               #ylab="mgp = c(3, 2, 0)",
               # font.label = list(size=10,face="bold",color='black',position="dodge",
               #                   vjust=0.5),
               #  y.label=0.7,
               # palette = palette,
               # 按照cyl填充颜色
               # palette = palette, # 修改颜色
               sorting = "descending",
               
               add = "segments",                             # 添加棒子
               rotate = TRUE,
               ggtheme = theme_pubr(),                        # 改变主题
               ylab="Correlation",
               xlab='',
               # dot.size = 6 ,
               title="RiskScore_ImmuneFunction_Correlation"
)+scale_colour_gradient( high = "yellow",low = "#246BAE")
p10<-p0+
  theme(legend.position = "right",
        panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
  theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
        axis.text.x =element_text(size=12,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=12,family = "Times", face = "bold"),
        plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
        legend.text = element_text(size = 16, family = "Times"),
        legend.title = element_text(size = 18, family = "Times",face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p10
ggsave("12_RiskScore_ImmuneFunction_Correlation.pdf",height = 6, width = 10, p10)
ggsave("12_RiskScore_ImmuneFunction_Correlation.png",height = 6, width = 10, p10)


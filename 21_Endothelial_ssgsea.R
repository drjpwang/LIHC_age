rm(list = ls())
setwd("/data/LIHC/")
if(! dir.exists("./21_Endothelial_ssgsea")){
  dir.create("./21_Endothelial_ssgsea")
}
setwd("./21_Endothelial_ssgsea")

getwd()


risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)
risk$sample <- rownames(risk)
risk1 <- risk[,c("sample","risk")]

#输入疾病样本表达矩阵fpkm
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)

cnames <- risk1$sample
train_t <- data[,cnames]

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
## 03-1 ssGSEA--------
library(GSVA)
library(Hmisc)
library(lance)
library(dplyr)
library(WGCNA)
library(ComplexHeatmap)
library(corrplot)
library(Hmisc)

load("../05_scRNA-seq/LIHC_sc_anno.rda")

cell_marker <- FindAllMarkers(sc_dat2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cell_Htop20 <- cell_marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少

colnames(cell_Htop20)
gene_set <- cell_Htop20[,c("cluster","gene")]

dat.final2 <- as.matrix(train_t)
gene_list <- split(as.matrix(gene_set)[,2],
                   gene_set[,1])

ssgsea_score = gsva(dat.final2, gene_list, 
                    method = "zscore", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.csv(ssgsea_score,
          file = "ssgsea_result.csv",quote = F)

# save.image("result.RData")
# load("result.RData")

mypalette <- colorRampPalette(brewer.pal(12,"Paired"))

risk <- risk[order(risk$risk),]

pdf('01.heatmap.pdf',  w=10,h=10,family='Times')
densityHeatmap(ssgsea_score ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = risk$risk, col = list(Group = c("High" = "red", "Low" = "blue"))) %v%
  Heatmap(ssgsea_score,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression",
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(10, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("blue", "white","red"))(100))
dev.off()

png('01.heatmap.png',w=10,h=10,units='in',res=600,family='Times')
densityHeatmap(ssgsea_score ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = risk$risk, col = list(Group = c("High" = "red", "Low" = "blue"))) %v%
  Heatmap(ssgsea_score,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression",
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(10, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("blue", "white","red"))(100))
dev.off()
## 堆积图-----------------
cibersort.res <- t(ssgsea_score)%>%as.data.frame()
res1 <- t(cibersort.res) %>% as.data.frame()
res1$cell_type <- rownames(res1)

# pdf('02.cell_barplotData.pdf',w=10,h=8, family = 'Times')
# res1 %>%
#   gather(id, fraction, -cell_type) %>%
#   merge(group,by='id')%>%
#   # 绘制堆积条形图
#   ggplot(aes(x=id, y=fraction, fill=cell_type)) +
#   geom_bar(position = 'stack',stat = 'identity')+
#   scale_y_continuous(expand = c(0,0))+
#   theme_bw()+
#   labs(x='',
#        y='Relative Percent',
#        fill='')+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = 'top') +
#   # scale_fill_manual(values = mypalette(28))+
#   facet_grid(~risk,scales= "free",space= "free")
# dev.off()
# 
# png('02.cell_barplotData.png',w=10,h=8,units='in',res=600, family = 'Times')
# res1 %>%
#   gather(id, fraction, -cell_type) %>%
#   merge(group,by='id')%>%
#   # 绘制堆积条形图
#   ggplot(aes(x=id, y=fraction, fill=cell_type)) +
#   geom_bar(position = 'stack',stat = 'identity')+
#   scale_y_continuous(expand = c(0,0))+
#   theme_bw()+
#   labs(x='',
#        y='Relative Percent',
#        fill='')+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = 'top') +
#   # scale_fill_manual(values = mypalette(28))+
#   facet_grid(~risk,scales= "free",space= "free")
# dev.off()


## 箱线图-------------
##  R值检验
### 过滤掉表达为0的
res2<-cibersort.res[,which(colSums(cibersort.res) > 0)] # 28
res2<-res2[rownames(risk),order(colnames(res2),decreasing = F)]
# group<-data.frame(sample=rownames(risk),Risk=risk$Risk)
identical(rownames(risk),rownames(res2))

### 按照字母进行排序
res2 <- res2 %>% as.data.frame %>% tibble::rownames_to_column(var = "id")
identical(res2$id,rownames(risk))
risk1 <- risk
risk1$id <- rownames(risk1)
risk1 <- risk1[,colnames(risk1) %in% c("id","risk")]
res2 <- merge(res2, risk1, by = "id")

### 删除OS、OStime、各基因和riskscore
# res2 <- res2[,-c(23:(length(colnames(res2))-2),length(colnames(res2)))]
colnames(res2)
## 差异检验
dat.test <- tidyr::gather(res2, Cell, `Cell composition`, -c(id, risk))
colnames(dat.test)
stat_res <- dat.test %>%
  group_by(Cell) %>%
  wilcox_test(`Cell composition` ~ risk) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_res
write.csv(stat_res, file = "03.cell_wilcoxon_res.csv", row.names = F, quote = F)

## 差异细胞
dif_cell<-stat_res[stat_res$p<0.05,]
dim(dif_cell)
dif_cell$Cell
# #过滤低表达细胞
# cibersort.res.fil <- data.frame(t(cibersort.res))
# m.mad <- apply(cibersort.res.fil,1,mad)
# cibersort.res.fil <- cibersort.res.fil[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2], 0.001)),]
# #既差异又高表达
# sigCell_use <- intersect(rownames(cibersort.res.fil),dif_cell$Cell)
# sigCell_use
#
# #保留sigCell_use数据，准备画图
# dat.test2 <- dat.test[dat.test$Cell == c('Macrophages M0','T cells CD8','T cells follicular helper'),]
# unique(dat.test2$Cell)
dat.test2 <- dat.test
dif_cell2 <- stat_res

## 绘图 ------------------------------------------------------------------------
#dif_cell2 <- dif_cell[dif_cell$Cell %in% sigCell_use,]
library(ggpubr)

p <- ggboxplot(dat.test2 , x = "Cell", y = "Cell composition",
               fill= "risk",
               palette = c("#dc143c", "#6495ed"),
               #outlier.shape = 0,
               outlier.size = 0.8,
               bxp.errorbar = T,
               #add = "mean_se"
               
) +ylim(0,1)+
  stat_pvalue_manual(dif_cell2 ,
                     x = "Cell",
                     y.position = 1,
                     size =4.5,
                     color = "black",
                     family = "Times",
                     # label = paste0("~italic(p)=={p.adj2}"),
                     label = "p.signif"
                     #parse = T
  ) +
  theme_bw()+
  labs(x = "", y = "Cell composition", color = "") +
  theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times"),
        axis.title.y = element_text(size = 20, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 13, color='black',face = "bold", family = "Times", angle = 60, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 15, face = "bold", family = "Times"),
        legend.text = element_text(size = 15, family = "Times",face = "bold"),
        legend.title = element_text(size = 18, family = "Times",face = "bold"),
        #plot.margin = ggplot2::margin(t=.3,b=0,l=2,r=.5, unit = "cm"),
        text = element_text(family = "Times"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.position = "top")
p
ggsave(filename = "04.Cell_boxplot.pdf", width = 10, height = 8)
ggsave(filename = "04.Cell_boxplot.png", width = 10, height = 8)



##相关性-----------------
res3 <- res2 %>% column_to_rownames(var = 'id')
num <- ncol(res3)
res3 <- res3[,-num]
name <- subset(dif_cell2,dif_cell2$p <0.05)
name <- name$Cell
res3 <- res3[,colnames(res3)%in%name]
stat.test <- corr.test(res3,method = 'spearman')
m <- stat.test$r
m <- t(m)
p <- stat.test$p
p<-t(p)
cor1 <- reshape2::melt(m)
p1 <- reshape2::melt(p)
colnames(p1)[3] <- 'pvalue'
both <- merge(cor1,p1,by = c('Var1','Var2'))
write.csv(both,file = '03.Correlaion.csv',row.names = T,quote = F)

tmp<- both

plot_cor<- m[rownames(m) %in% unique(tmp$Var1),]
plot_cor<- plot_cor[,unique(tmp$Var2)]
plot_p<- p[rownames(p) %in% unique(tmp$Var1),]
plot_p<- plot_p[,unique(tmp$Var2)]

pdf('03.Correlation.pdf',width = 18,height = 12)
par(mar = c(5, 10, 2, 2))
corrplot(plot_cor, p.mat =  plot_p,type = 'upper',
         insig = 'label_sig',
         sig.level = c(0.001, 0.01,0.05),
         col=rev(COL2('RdYlBu', 200)),tl.pos="tp",
         #sig.level = -1,
         pch.cex = 1, pch.col = 'grey20')
corrplot(corr = plot_cor, 
         add=TRUE, type="lower", method="number",
         col='black',diag=FALSE,
         tl.pos="n", cl.pos="n",number.cex = 1)

dev.off()
png('03.Correlation.png',width = 1800,height = 1000)
par(mar = c(5, 10, 2, 2))
corrplot(plot_cor, p.mat =  plot_p,type = 'upper',
         insig = 'label_sig',
         sig.level = c(0.001, 0.01,0.05),
         col=rev(COL2('RdYlBu', 200)),tl.pos="tp",
         #sig.level = -1,
         pch.cex = 1, pch.col = 'grey20')
corrplot(corr = plot_cor, 
         add=TRUE, type="lower", method="number",
         col='black',diag=FALSE,
         tl.pos="n", cl.pos="n",number.cex = 1)

dev.off()

##候选基因-------------------
names <- stat_res[stat_res$p < 0.05,]
names <- names$Cell
diff_score <- t(res3[,names])
hubgene <- readRDS("../05_Survival/08_lasso_geneids.rds")
hubgene <- as.data.frame(hubgene)
colnames(hubgene) <- 'symbol'
# hubgene <- hubgene[-3,]%>%as.data.frame()
hub_exp <- train_t[hubgene$symbol,]
hub_exp <- t(hub_exp) %>% as.data.frame()

diff_score <- t(diff_score) %>% as.data.frame()

#比较两个对象是否相等
identical(rownames(hub_exp), rownames(diff_score))
# expr <- expr[rownames(expr) %in% rownames(tiics_result),]
expr <- hub_exp
tiics_result <- diff_score

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
write.csv(cor_dat,"04_DECell&Gene_correlation_cor.csv")

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
ggsave("05_cibersort_Biomarker_ImmuneCell_corplot.pdf", ciber_corplot3, width = 10, height = 6)
ggsave("05_cibersort_Biomarker_ImmuneCell_corplot.png", ciber_corplot3, width = 10, height = 6)


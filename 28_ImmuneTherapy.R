rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./28_ImmuneTherapy")){
  dir.create("./28_ImmuneTherapy")
}
setwd("./28_ImmuneTherapy")

getwd()


library(GEOquery)

#下载训练集数据
geo_id1 <- "GSE140901"
expr <- read.table("/data/home/LIHC/GEO_immune/GSE140901/GSE140901.txt", row.names = 1, header = T)
range(expr)

group_t <- read.table("/data/home/LIHC/GEO_immune/GSE140901/group.txt", header = T)
table(group_t$group)

identical(group_t$sample, colnames(expr))

train_expr <- expr

write.csv(train_expr, paste0("01_",geo_id1,"_train_expr.csv"))
save(train_expr, file=paste0("01_",geo_id1,"_train_expr.rda"))

rownames(group_t) <- group_t$sample

write.csv(group_t, paste0("02_",geo_id1,"_train_infor.csv"))
save(group_t, file=paste0("02_",geo_id1,"_train_infor.rda"))


rm(list = ls())
geo_id1 <- "GSE181946"
expr <- read.table("/data/home/LIHC/GEO_immune/GSE181946/GSE181946_tpm.txt", row.names = 1, header = T)
range(expr)

group_t <- read.table("/data/home/LIHC/GEO_immune/GSE181946/group.txt", header = T)
table(group_t$group)

identical(group_t$sample, colnames(expr))

expr <- expr[,group_t$sample]

train_expr <- expr

write.csv(train_expr, paste0("01_",geo_id1,"_train_expr.csv"))
save(train_expr, file=paste0("01_",geo_id1,"_train_expr.rda"))

rownames(group_t) <- group_t$sample

write.csv(group_t, paste0("02_",geo_id1,"_train_infor.csv"))
save(group_t, file=paste0("02_",geo_id1,"_train_infor.rda"))



rm(list = ls())
geo_id1 <- "GSE202069"
expr <- read.table("/data/home/LIHC/GEO_immune/GSE202069/GSE202069_gene_tpm_expression.txt", row.names = 1, header = T)
range(expr)

group_t <- read.table("/data/home/LIHC/GEO_immune/GSE202069/group.txt", header = T)
table(group_t$group)

identical(group_t$sample, colnames(expr))

expr <- expr[,group_t$sample]

train_expr <- expr

write.csv(train_expr, paste0("01_",geo_id1,"_train_expr.csv"))
save(train_expr, file=paste0("01_",geo_id1,"_train_expr.rda"))

rownames(group_t) <- group_t$sample

write.csv(group_t, paste0("02_",geo_id1,"_train_infor.csv"))
save(group_t, file=paste0("02_",geo_id1,"_train_infor.rda"))




rm(list = ls())
geo_id1 <- "GSE215011"
expr <- read.table("/data/home/LIHC/GEO_immune/GSE215011/GSE215011_fpkm.txt", row.names = 1, header = T)
range(expr)

group_t <- read.table("/data/home/LIHC/GEO_immune/GSE215011/group.txt", header = T)
table(group_t$group)

identical(group_t$sample, colnames(expr))

expr <- expr[,group_t$sample]
expr <- log2(expr+1)

train_expr <- expr

write.csv(train_expr, paste0("01_",geo_id1,"_train_expr.csv"))
save(train_expr, file=paste0("01_",geo_id1,"_train_expr.rda"))

rownames(group_t) <- group_t$sample

write.csv(group_t, paste0("02_",geo_id1,"_train_infor.csv"))
save(group_t, file=paste0("02_",geo_id1,"_train_infor.rda"))


rm(list = ls())
if(! dir.exists("./Cell_ImmuneTherapy")){
  dir.create("./Cell_ImmuneTherapy")
}
setwd("./Cell_ImmuneTherapy")

getwd()

#计算响应差异
#关注细胞的通路富集
# sc_dat2 <- readRDS("../../25_T_Endo/ALL_celltype_scdata.RDS")
load("../../25_Tcell/T_sc_anno.rda")

cell.markers <- FindAllMarkers(T_dat2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cell.markers, "01_cell.markers.csv", quote = F, row.names = F)


cell_Htop100 <- cell.markers %>% group_by(cluster) %>% top_n(n = 150, wt = avg_log2FC) ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少
Endo <- cell_Htop100[cell_Htop100$cluster == "Endothelial",]
Treg <- cell_Htop100[cell_Htop100$cluster == "Treg",]


marker1 <- Endo$gene
marker1 <- c(marker1,"PECAM1","VWF","CLDN5","RAMP2")
marker1 <- unique(marker1)
cell1 <- c(rep("Endothelial",length(marker1)))

marker_gene1 <- data.frame(marker=marker1, celltype=cell1)


marker2 <- Treg$gene
cell2 <- c(rep("Treg",nrow(Treg)))

marker_gene2 <- data.frame(marker=marker2, celltype=cell2)

marker_gene <- rbind(marker_gene1, marker_gene2)


load("../01_GSE202069_train_expr.rda")

load("../02_GSE202069_train_infor.rda")

train_t <- train_expr
#计算细胞的ssgsea评分
library(GSVA)
gene_list <- split(as.matrix(marker_gene)[,1], marker_gene[,2])
ssgsea_score1 <- gsva(as.matrix(train_t),  gene_list,
                      method = "zscore", kcdf='Gaussian',abs.ranking=TRUE) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

write.csv(ssgsea_score1, "02_Cell_ssgsea_score.csv")


#患者响应与不响应样本样本之间的细胞浸润差异箱式图
sig_res2 <- ssgsea_score1
sig_cell <- as.data.frame(row.names(sig_res2))

colnames(sig_cell) <- "Cell"

sig_res2 <- cbind(sig_cell, sig_res2)
rownames(sig_res2) <- 1:dim(sig_res2)[1]
write.csv(sig_res2, "03_cell_result_cell.csv")

# keep <- rowSums(sig_res2>0)>=floor(0.3*ncol(sig_res2))  # 过滤掉在所有样本中占比低于30%的细胞类型
# sig_res2 <- sig_res2[keep,]


dat.cibersort <- sig_res2 %>% 
  #tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
colnames(dat.cibersort) <- dat.cibersort[1,]
dat.cibersort <- dat.cibersort[-1,]
colnames(dat.cibersort)[1] <- "sample"

dat.cibersort <- merge(group_t, dat.cibersort, by = "sample")

# dat.cibersort <- subset(dat.cibersort, select=-c(2))
dat.cibersort1 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))

colnames(dat.cibersort1)
library(rstatix)

dat.cibersort1$Score <- as.numeric(dat.cibersort1$Score)
stat_cibersort <- dat.cibersort1 %>% 
  group_by(ImmuneCell) %>%
  rstatix::wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.csv(stat_cibersort,file = '04_LIHC_Cell.csv',row.names = T)

DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort,file = '05_ImmuneTherapy_DE.csv',row.names = T)

colnames(dat.cibersort1)
#violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
violin.cibersort<-dat.cibersort1
library(ggpubr)
ciber_boxplot <- ggplot(dat.cibersort1, aes(x=ImmuneCell,
                                            y=Score,
                                            fill=group)) +
  # geom_violin(trim=F,color="black") + "#FF7F24","#5F9EA0"
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values= c("red","grey"), name = "Group")+
  labs(title="Immune Therapy", x="", y = "Fraction",size=20) +
  stat_compare_means(data = dat.cibersort1,
                     mapping = aes(group = group),
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
ggsave("06_cellImmuneTherapy_boxplot.pdf", ciber_boxplot, width = 8, height = 6)
ggsave("06_cellImmuneTherapy_boxplot.png", ciber_boxplot, width = 8, height = 6)


data <- t(ssgsea_score1) %>% as.data.frame()
data1 <- data

p1<-ggplot(data1,aes(x=Endothelial,y=Treg))+geom_point(aes(x=Endothelial,y=Treg), size=2,alpha=0.9)+
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
  labs(x="Endothelail",y="Treg")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1

pdf('06_Endothelail_Treg_Correlation.pdf',width=7,height=6,family='Times')
print(p1)
dev.off()

png('06_Endothelail_Treg_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p1)
dev.off()# stat.test <- corr.test(res3,method = 'spearman')

ssgsea_score2 <- ssgsea_score1


#内皮细胞亚型
load("../../21_Endothelial/Endothelial_sc_anno.rda")

endo.markers <- FindAllMarkers(sc_dat2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(endo.markers, "01_endo.markers.csv", quote = F, row.names = F)

# a <- "IGFBP3"

endo_Htop100 <- endo.markers %>% group_by(cluster) %>% top_n(n = 150, wt = avg_log2FC) ##n=10 根据avg_logFC这一列把前10个差异基因取出来，小于的话，有多少取多少
# AE1 <- endo_Htop100[endo_Htop100$cluster == "Artery-1 ECs",]
# setdiff(a, AE1$gene)  # 排名18

# marker <- endo_Htop100$gene
# cell <- c(rep("Artery-1 ECs",nrow(AE1)))

marker_gene <- data.frame(marker=endo_Htop100$gene, celltype=endo_Htop100$cluster)


load("../01_GSE202069_train_expr.rda")

load("../02_GSE202069_train_infor.rda")

train_t <- train_expr
#计算细胞的ssgsea评分
library(GSVA)
gene_list <- split(as.matrix(marker_gene)[,1], marker_gene[,2])
ssgsea_score1 <- gsva(as.matrix(train_t),  gene_list,
                      method = "zscore", kcdf='Gaussian',abs.ranking=TRUE) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

write.csv(ssgsea_score1, "07_Cell_ssgsea_score.csv")


#患者响应与不响应样本样本之间的细胞浸润差异箱式图
sig_res2 <- ssgsea_score1
sig_cell <- as.data.frame(row.names(sig_res2))

colnames(sig_cell) <- "Cell"

sig_res2 <- cbind(sig_cell, sig_res2)
rownames(sig_res2) <- 1:dim(sig_res2)[1]
write.csv(sig_res2, "08_Endothelial_result_cell.csv")

# keep <- rowSums(sig_res2>0)>=floor(0.3*ncol(sig_res2))  # 过滤掉在所有样本中占比低于30%的细胞类型
# sig_res2 <- sig_res2[keep,]


dat.cibersort <- sig_res2 %>% 
  #tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
colnames(dat.cibersort) <- dat.cibersort[1,]
dat.cibersort <- dat.cibersort[-1,]
colnames(dat.cibersort)[1] <- "sample"

dat.cibersort <- merge(group_t, dat.cibersort, by = "sample")

# dat.cibersort <- subset(dat.cibersort, select=-c(2))
dat.cibersort1 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))

colnames(dat.cibersort1)
library(rstatix)

dat.cibersort1$Score <- as.numeric(dat.cibersort1$Score)
stat_cibersort <- dat.cibersort1 %>% 
  group_by(ImmuneCell) %>%
  rstatix::wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.csv(stat_cibersort,file = '09_LIHC_Endothelial.csv',row.names = T)

DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort,file = '10_Endothelial_ImmuneTherapy_DE.csv',row.names = T)

colnames(dat.cibersort1)
#violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
violin.cibersort<-dat.cibersort1
library(ggpubr)
ciber_boxplot <- ggplot(dat.cibersort1, aes(x=ImmuneCell,
                                            y=Score,
                                            fill=group)) +
  # geom_violin(trim=F,color="black") + "#FF7F24","#5F9EA0"
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values= c("red","grey"), name = "Group")+
  labs(title="Immune Therapy", x="", y = "Fraction",size=20) +
  stat_compare_means(data = dat.cibersort1,
                     mapping = aes(group = group),
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
ggsave("11_cellImmuneTherapy_boxplot.pdf", ciber_boxplot, width = 8, height = 6)
ggsave("11_cellImmuneTherapy_boxplot.png", ciber_boxplot, width = 8, height = 6)


data <- t(ssgsea_score1) %>% as.data.frame()
ssgsea_score2 <- t(ssgsea_score2) %>% as.data.frame()
identical(rownames(data), rownames(ssgsea_score2))
colnames(data)
data1 <- cbind(data, ssgsea_score2)
data2 <- data1[,c("Artery-1 ECs","Treg")]

p2<-ggplot(data2,aes(x=`Artery-1 ECs`,y=Treg))+geom_point(aes(x=`Artery-1 ECs`,y=Treg), size=2,alpha=0.9)+
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
  labs(x="Artery-1 ECs",y="Treg")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2

pdf('11_Artery-1 ECs_Treg_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('11_Artery-1 ECs_Treg_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()# stat.test <- corr.test(res3,method = 'spearman')



setwd("../")

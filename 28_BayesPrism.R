rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./28_BayesPrism")){
  dir.create("./28_BayesPrism")
}
setwd("./28_BayesPrism")

getwd()


if(! dir.exists("./GSE181946")){
  dir.create("./GSE181946")
}
setwd("./GSE181946")

library(SingleR)
library(survminer)
library(survival)
library(Seurat)

library("devtools")
install_github("Danko-Lab/BayesPrism/BayesPrism")
library(BayesPrism)

load("../../21_Endothelial/Endothelial_sc_anno.rda")

cell.type.labels = sc_dat2$cell_type
cell.state.labels = NULL
cell_exp <- as.matrix(GetAssayData(sc_dat2, layer = "count", assay = "RNA"))


#TCGA数据
bulk_exp = read.table('/data/home/LIHC/GEO_immune/GSE181946/GSE181946_count.txt',
                      row.names = 1,header = T)
range(bulk_exp)
bulk_exp <- round(bulk_exp)
colnames(bulk_exp) <- gsub("\\.","\\-",colnames(bulk_exp))

group <- read.csv('../../28_ImmuneTherapy/02_GSE181946_train_infor.csv',header = T,row.names = 1)
# Tumor <- group[group$group%in%'GSC',]$sample
# bulk_exp <- bulk_exp[,Tumor]
colnames_order <- str_sort(unique(sc_dat2$cell_type), numeric = TRUE)
# common <- intersect(rownames(bulk_exp), rownames(cell_exp))
# drop   <- grepl("^MT-|^RPS|^RPL|^RP[SL]|^HB[AB]|^IG[HKL]|^LINC", common, ignore.case = TRUE)
# genes  <- setdiff(common, common[drop])
# cell_exp  <- cell_exp[genes, , drop = FALSE]

myPrism <- new.prism(
  reference = t(cell_exp),  ##行是细胞
  mixture = t(bulk_exp),##行是样本
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL  # cell.type.labels中恶性细胞亚群
)


bp.res <- run.prism(
  prism = myPrism,
  n.cores = 100
)

save(bp.res, file = 'ImmuneTherapy_bp.res.RData')

# load('bp.res.RData')
theta_type<- get.fraction(bp=bp.res,
                          which.theta="final",##initial 开始的初始估计，all返回每次迭代的全部结果
                          state.or.type="type")##type 按细胞类型汇总的比例，state 按细胞状态或亚群
# theta_state <- bp.res@posterior.initial.cellState@theta # 只能在第一次运行结束的时候提取先验亚群的反卷积信息
write.csv(theta_type,file = 'ImmuneTherapy_Bayes.csv')


identical(rownames(theta_type), group$sample)
bp_result <- as.data.frame(theta_type)
bp_result$sample <- rownames(bp_result)

Bayes_result <- merge(bp_result, group, by = "sample")
Bayes_result1 <- tidyr::gather(Bayes_result, ImmuneCell, Score, -c("sample", "group"))
ciber_boxplot <- ggplot(Bayes_result1, aes(x=ImmuneCell,
                                           y=Score,
                                           fill=group)) +
  # geom_violin(trim=F,color="black") + "#FF7F24","#5F9EA0"
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values= c("#dc143c","#6495ed"), name = "Group")+
  labs(title="Immune Therapy", x="", y = "Fraction",size=20) +
  stat_compare_means(data = Bayes_result1,
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
ggsave("01_cellBayes_ImmuneTherapy_boxplot.pdf", ciber_boxplot, width = 8, height = 6)
ggsave("01_cellBayes_ImmuneTherapy_boxplot.png", ciber_boxplot, width = 8, height = 6)


#K-M曲线
data=read.table("../../00.data/TCGA/TCGA-LIHC-mRNA-counts.txt",sep="\t",header=T,row.names = 1,
                stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)
range(data)

load("../../03_Survival/train_riskScore.rda")

train_t <- data[, rownames(risk)]

myPrism <- new.prism(
  reference = t(cell_exp),  ##行是细胞
  mixture = t(train_t),##行是样本
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL  # cell.type.labels中恶性细胞亚群
)


bp.res <- run.prism(
  prism = myPrism,
  n.cores = 100
)

save(bp.res, file = 'KM_bp.res.RData')

# load('bp.res.RData')
theta_type<- get.fraction(bp=bp.res,
                          which.theta="final",##initial 开始的初始估计，all返回每次迭代的全部结果
                          state.or.type="type")##type 按细胞类型汇总的比例，state 按细胞状态或亚群
# theta_state <- bp.res@posterior.initial.cellState@theta # 只能在第一次运行结束的时候提取先验亚群的反卷积信息
write.csv(theta_type,file = 'KM_Bayes.csv')


clinical <- risk[,c("OS","OS.time")]

bp_score <- as.data.frame(theta_type)
bp_score$sample <- rownames(bp_score)
bp_score$sample <- gsub("\\.", "_", bp_score$sample)
bp_score$sample <- gsub("-", "_", bp_score$sample)

surv <- clinical[rownames(bp_score),]
surv$sample <- rownames(surv)
surv$sample <- gsub("-", "_", surv$sample)


#整合表达谱数据和生存信息
train_v <- merge(bp_score, surv, by = "sample")

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
#最佳截断值
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "Artery_1_ECs",  #需要计算的数据列名
                                   minprop  = 0.3
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res # 0.0007308237

diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"Artery_1_ECs"] > res["Artery_1_ECs",1],'High','Low')
table(diff_expr_clinical1$risk1)

# diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"Artery_1_ECs"] > median(diff_expr_clinical1[,"Artery_1_ECs"]),'High','Low')
# table(diff_expr_clinical1$risk1)

library(survminer)
library(survival)
# K-M曲线 ----------
range(diff_expr_clinical1$OS.time)
diff_expr_clinical1$OS.time <- diff_expr_clinical1$OS.time/30

kmfit <- survfit(Surv(OS.time, OS) ~ risk1, data = diff_expr_clinical1)
surv_pvalue(kmfit) # GB TW PP mPP FH
pmethod <- surv_pvalue(kmfit, method = "GB")
pmethod


png('02_KM_Artery_1_ECs.png',w=7,h=7,family='Times',units = 'in',res = 600)
s_surv <- ggsurvplot(kmfit,
                     pval = T,   #!!
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

pdf('02_KM_Artery_1_ECs.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
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



rm(list = ls())
load("../25_Tcell/T_sc_anno.rda")

cell.type.labels = T_dat2$cell_type
cell.state.labels = NULL
cell_exp <- as.matrix(GetAssayData(T_dat2, layer = "count", assay = "RNA"))

#TCGA数据
bulk_exp = read.table('/data/home/LIHC/GEO_immune/GSE181946/GSE181946_count.txt',
                      row.names = 1,header = T)
range(bulk_exp)
bulk_exp <- round(bulk_exp)
colnames(bulk_exp) <- gsub("\\.","\\-",colnames(bulk_exp))

group <- read.csv('../../28_ImmuneTherapy/02_GSE181946_train_infor.csv',header = T,row.names = 1)
# Tumor <- group[group$group%in%'GSC',]$sample
# bulk_exp <- bulk_exp[,Tumor]
colnames_order <- str_sort(unique(T_dat2$cell_type), numeric = TRUE)
# common <- intersect(rownames(bulk_exp), rownames(cell_exp))
# drop   <- grepl("^MT-|^RPS|^RPL|^RP[SL]|^HB[AB]|^IG[HKL]|^LINC", common, ignore.case = TRUE)
# genes  <- setdiff(common, common[drop])
# cell_exp  <- cell_exp[genes, , drop = FALSE]

myPrism <- new.prism(
  reference = t(cell_exp),  ##行是细胞
  mixture = t(bulk_exp),##行是样本
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL  # cell.type.labels中恶性细胞亚群
)


bp.res <- run.prism(
  prism = myPrism,
  n.cores = 100
)

save(bp.res, file = 'Allcell.bp.res.RData')

# load('bp.res.RData')
theta_type<- get.fraction(bp=bp.res,
                          which.theta="final",##initial 开始的初始估计，all返回每次迭代的全部结果
                          state.or.type="type")##type 按细胞类型汇总的比例，state 按细胞状态或亚群
# theta_state <- bp.res@posterior.initial.cellState@theta # 只能在第一次运行结束的时候提取先验亚群的反卷积信息
write.csv(theta_type,file = 'Allcell.Bayes.csv')


identical(rownames(theta_type), group$sample)
bp_result <- as.data.frame(theta_type)
bp_result$sample <- rownames(bp_result)

Bayes_result <- merge(bp_result, group, by = "sample")
Bayes_result1 <- tidyr::gather(Bayes_result, ImmuneCell, Score, -c("sample", "group"))
ciber_boxplot <- ggplot(Bayes_result1, aes(x=ImmuneCell,
                                           y=Score,
                                           fill=group)) +
  # geom_violin(trim=F,color="black") + "#FF7F24","#5F9EA0"
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values= c("#dc143c","#6495ed"), name = "Group")+
  labs(title="Immune Therapy", x="", y = "Fraction",size=20) +
  stat_compare_means(data = Bayes_result1,
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
ggsave("03_cellBayes_Allcell_boxplot.pdf", ciber_boxplot, width = 8, height = 6)
ggsave("03_cellBayes_Allcell_boxplot.png", ciber_boxplot, width = 8, height = 6)



#相关性分析
bp_data <- read.csv("ImmuneTherapy_Bayes.csv", header = T, row.names = 1)
bp_data$sample <- rownames(bp_data)

identical(Bayes_result$sample, bp_data$sample)
data1 <- merge(Bayes_result, bp_data, by = "sample")
colnames(data1)
colnames(data1) <- gsub("\\.", " ", colnames(data1))
colnames(data1)[c(21,23)] <- c("Artery-2 ECs","Artery-1 ECs")
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

pdf('04_Artery-1 ECs_Treg_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('04_Artery-1 ECs_Treg_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()



rm(list = ls())
sc_dat2 <- readRDS("../../25_T_Endo/ALL_celltype_scdata.RDS")

cell.type.labels = sc_dat2$subcluster
cell.state.labels = NULL
cell_exp <- as.matrix(GetAssayData(sc_dat2, layer = "count", assay = "RNA"))

#TCGA数据
data=read.table("../../00.data/TCGA/TCGA-LIHC-mRNA-counts.txt",sep="\t",header=T,row.names = 1,
                stringsAsFactors = F,quote="",check.names = F)
range(data)
data <- round(data)
colnames(data)=substring(colnames(data),1,12)

load("../../03_Survival/11_all_train_riskScore.rda")

train_t <- data[, rownames(risk)]


myPrism <- new.prism(
  reference = t(cell_exp),  ##行是细胞
  mixture = t(train_t),##行是样本
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL  # cell.type.labels中恶性细胞亚群
)


bp.res <- run.prism(
  prism = myPrism,
  n.cores = 100
)

save(bp.res, file = 'Allcell.Integration.bp.res.RData')

# load('bp.res.RData')
theta_type<- get.fraction(bp=bp.res,
                          which.theta="final",##initial 开始的初始估计，all返回每次迭代的全部结果
                          state.or.type="type")##type 按细胞类型汇总的比例，state 按细胞状态或亚群
# theta_state <- bp.res@posterior.initial.cellState@theta # 只能在第一次运行结束的时候提取先验亚群的反卷积信息
write.csv(theta_type,file = 'Allcell.Integration.Bayes.csv')


bp_result <- as.data.frame(theta_type)
# bp_result$sample <- rownames(bp_result)


#相关性分析
data1 <- bp_result
colnames(data1)
# colnames(data1) <- gsub("\\.", " ", colnames(data1))
# colnames(data1)[c(5,14)] <- c("Artery-2 ECs","Artery-1 ECs")
data2 <- data1[,c("Artery-1 ECs","Treg","Macrophage")]

p2<-ggplot(data2,aes(x=`Artery-1 ECs`,y=Treg))+geom_point(aes(x=`Artery-1 ECs`,y=Treg), size=2,alpha=0.9)+
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
  labs(x="Artery-1 ECs",y="Treg")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2

pdf('05_Artery-1 ECs_Macrophage_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('05_Artery-1 ECs_Macrophage_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()

###剔除离群样本
data2_1 <- data2[data2$Treg < 0.015,]
p2<-ggplot(data2_1,aes(x=`Artery-1 ECs`,y=Treg))+geom_point(aes(x=`Artery-1 ECs`,y=Treg), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=data2_1, method = "spearman",digits = 2)+
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

pdf('05_Artery-1 ECs_Treg_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('05_Artery-1 ECs_Treg_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()

hubgene <- c("IGFBP3","ELF1","TCF7L1","SPP1")

df=read.table("../../00.data/TCGA/TCGA-LIHC.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(4,3,2)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]

data=read.table("../../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[hubgene,]

colnames(data)=substring(colnames(data),1,12)
data=t(data)
sameSample=intersect(row.names(data),row.names(clinical))##365
data=data[sameSample,]


data3 <- as.data.frame(data)

identical(rownames(data3), rownames(data2))

##热图---
library(psych)
library(reshape2)
r_result <- corr.test(data3,data2,method = 'spearman')

r <-
  r_result$r %>%
  reshape2::melt() %>%
  set_names(c('tax', 'Index', 'r'))

p <-
  r_result$p %>%
  reshape2::melt() %>%
  set_names(c('tax', 'Index', 'P_value')) %>%
  mutate(P_value_sig = case_when(P_value > 0.05 ~ "ns",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))


data <- cbind(r,p) %>% dplyr::select(-(4:5)) #删除重复列
head(data)
write.csv(data,"06_cor_cell_gene.csv")
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

immune_corplot <- ggplot(data,aes(x=tax, y=Index))+
  geom_tile(aes(fill=r), color = 'white',alpha = 0.6) +
  geom_text(label = textMatrix, color = 'black', size =6) +
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
        axis.text.x = element_text(hjust = 1, size = 14, face = "bold", angle=45),
        axis.text.y = element_text(size = 14, face = "bold")) + 
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) + 
  # scale_x_discrete(position = "top") + 
  xlab('') +
  ylab('')
# ggsave('heatmap_corr.png',height = 2, width = 2)
immune_corplot
ggsave(filename = '07_Biomarker_ImmuneCell_corplot.pdf',immune_corplot,w=8,h=6)
ggsave(filename = '07_Biomarker_ImmuneCell_corplot.png',immune_corplot,w=8,h=6)

data2$sample <- rownames(data2)
data3$sample <- rownames(data3)
cell_gene <- merge(data2,data3,by="sample")
data4 <- cell_gene[,c("Artery-1 ECs","ELF1")]
data4 <- as.data.frame(data4)
rownames(data4) <- rownames(cell_gene)

p2<-ggplot(data4,aes(x=`Artery-1 ECs`,y=ELF1))+geom_point(aes(x=`Artery-1 ECs`,y=ELF1), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=data4, method = "spearman",digits = 2)+
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
  labs(x="Artery-1 ECs",y="ELF1")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2

pdf('08_Artery-1 ECs_ELF1_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('08_Artery-1 ECs_ELF1_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()



data5 <- cell_gene[,c("Artery-1 ECs","SPP1")]
data5 <- as.data.frame(data5)
rownames(data5) <- rownames(cell_gene)

p2<-ggplot(data5,aes(x=`Artery-1 ECs`,y=SPP1))+geom_point(aes(x=`Artery-1 ECs`,y=SPP1), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=data5, method = "spearman",digits = 2)+
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
  labs(x="Artery-1 ECs",y="SPP1")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2

pdf('09_Artery-1 ECs_SPP1_Correlation.pdf',width=7,height=6,family='Times')
print(p2)
dev.off()

png('09_Artery-1 ECs_SPP1_Correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p2)
dev.off()



data6 <- cell_gene[,c("IGFBP3","ELF1")]
data6 <- as.data.frame(data6)
rownames(data6) <- rownames(cell_gene)

p2<-ggplot(data6,aes(x=`IGFBP3`,y=ELF1))+geom_point(aes(x=`IGFBP3`,y=ELF1), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=data6, method = "spearman",digits = 2)+
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
p2

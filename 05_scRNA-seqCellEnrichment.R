rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./05_scRNA-seq_CellEnrichment")){
  dir.create("./05_scRNA-seqCellEnrichment")
}
setwd("./05_scRNA-seqCellEnrichment")

getwd()

load("../05_scRNA-seq/LIHC_sc_anno.rda")

pacman::p_load(reshape2,tidyverse,Seurat,rstatix,ggpubr,gridExtra)
scRNA.norm.umap <- sc_dat2
table(scRNA.norm.umap$cell_type)

# load("../02_CandidateGene/candidate_gene.rda")
load("../03_Survival/LASSO_GeneID.rda")

sce.all =  AddModuleScore(object = scRNA.norm.umap,list(lasso_geneids$symbol))
head(sce.all@meta.data)
stat.test <- sce.all@meta.data%>%
  group_by(cell_type)%>%
  wilcox_test(Cluster1 ~ type)%>%
  adjust_pvalue(method = 'fdr')

stat.test$p <- ifelse(stat.test$p<0.001,"***",
                      ifelse(stat.test$p<0.01,"**",
                             ifelse(stat.test$p<0.05,"*",'ns')))
temp1 <- stat.test[stat.test$p < 0.05, ]
write.csv(temp1,"temp1.csv")
table(sce.all@meta.data$sample_type)

External_va <- ggplot(sce.all@meta.data, aes(x = cell_type, y = Cluster1, color = type)) +
  # ylim(0,0.8)+
  stat_boxplot(geom="errorbar",
               width=0.2,
               position = position_dodge(0.9)) +
  geom_boxplot(aes(color=type,fill=type),
               width=0.8,
               position=position_dodge(0.9),
               outlier.shape = NA  )+
  scale_fill_manual(values= c("#4682B4","#CD3700"))+#绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  #
  # annotate("rect", xmin = 0.4, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#c1f1fc") +
  #  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#ebffac") +
  #添加散点
  # geom_dotplot(dotsize = 0.4,binaxis = "y", stackdir = "center",position = position_dodge(0.8))+
  #添加分组辅助线
  #geom_vline(xintercept = 1.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  #geom_vline(xintercept = 2.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  #geom_vline(xintercept = 3.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  #geom_vline(xintercept = 4.5, lty="dashed", color = "grey50", linewidth = 0.8)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12),
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title=" ", x="", y = "AddModuleScore",size=20) +
  scale_colour_manual(values = c("black","black"))+
  #   scale_color_brewer(palette = 'Set2')+
  stat_compare_means(aes(group = type,label = "p.format"),
                     method = "wilcox.test",label = 'p.signif')#+,label.x = 1.4
# stat_pvalue_manual(stat.test,
#             y.position = 0.2,
#            size = 3.2,
#            family = "Times",
#            label = "p",
#            parse = T,
#            face = "bold")

External_va
ggsave(filename = '01_AddModuleScore_diff.pdf', External_va, w=8, h=6)
ggsave(filename = '01_AddModuleScore_diff.png', External_va, w=8, h=6, dpi = 600)



train_plot <- ggplot(sce.all@meta.data, aes(x=cell_type, y=Cluster1, fill = type)) +
  scale_fill_manual(values = c("#FF7F24","#5F9EA0")) + 
  #ylim(4.2,11.3)+
  geom_boxplot()+
  stat_compare_means(aes(group = type),method = 'wilcox.test',label = "p.signif",cex=4)+
  # annotate(geom = "text", x = violin_dat1$gene, y = 15, size = 5,
  # label =as.character( violin_dat1$P.Value.signif)) +
  xlab("Cell Type") + 
  ylab("AddModuleScore") + 
  # labs(title = "AddModuleScore")+ 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+
  theme(axis.text.x=element_text(face = "bold", color="black",angle = 70,vjust = 1, hjust = 1 ))+#angle = 70
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size=16),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size=12))+
  theme(plot.title = element_text(size=16))
train_plot
ggsave(filename = '13.AddModuleScore_diff1.pdf', train_plot, w=8, h=6)
ggsave(filename = '13.AddModuleScore_diff1.png', train_plot, w=8, h=6, dpi = 600)


marker <- c(#"EPCAM","KRT19","CDH1", # epithelial "KRT18",
  "DCN","THY1","COL1A1","COL1A2", # fibroblast
  "PECAM1","VWF","CLDN5","RAMP2", # endothelial
  # "ACTA2","TAGLN",         # Smooth Muscle cell
  # "ALB","APOA2","TTR",            # Hepatocytes
  # "RGS5","MYL9","KCNJ8",           # Pericytes
  # "ITGAX",                         # Myeloid
  "CD3D","CD3E","CD3G","TRAC",    # T cell
  "NKG7","GNLY","KLRD1",        # NK cell
  "CD3D","CD3E","CD3G","TRAC",    # T cell
  "NKG7","GNLY","KLRD1",        # NK cell
  "IGHM","IGHG3","IGHA2","CD79A", # B cell
  "CD79A","MS4A1","CD19",      # Plasma
  "LYZ","CLEC9A","FCER1A",        # DC  "CCR7",
  # "CSF3R","FCGR3B",               # Neutrophil
  # "KIT","MS4A2","GATA2",          # Mast cell
  "CD163","CCL18","CD68")  # macrophage  ,"CD86"
cell <- c("Fibroblast","Fibroblast","Fibroblast","Fibroblast",
          "Endothelial","Endothelial","Endothelial","Endothelial",
          "T cell","T cell","T cell","T cell",
          "NK cell","NK cell","NK cell",
          "NKT cell","NKT cell","NKT cell","NKT cell","NKT cell","NKT cell","NKT cell",
          "B cell","B cell","B cell","B cell",
          "Plasma","Plasma","Plasma",
          "DC","DC","DC",
          "Macrophage","Macrophage","Macrophage")
marker_gene <- data.frame(marker=marker, celltype=cell)


data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)

load("../03_Survival/train_riskScore.rda")

train_t <- data[, rownames(risk)]
gene_list <- split(as.matrix(marker_gene)[,1], marker_gene[,2])


#计算细胞的ssgsea评分
library(GSVA)
ssgsea_score1 <- gsva(as.matrix(train_t),  gene_list,
                      method = "ssgsea", kcdf='Gaussian',abs.ranking=TRUE) # method = "ssgsea"，"gsva"，"zscore"，"plage"无显著差异
# "Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts

write.csv(ssgsea_score1, "02_Cell_ssgsea_score.csv")

#读取生存信息
df=read.table("../00.data/TCGA/TCGA-LIHC.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(4,3,2)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]


ssgsea_score1 <- t(ssgsea_score1) %>% as.data.frame()
surv <- clinical[rownames(ssgsea_score1),]


ssgsea_score2 <- ssgsea_score1
ssgsea_score2$sample <- rownames(ssgsea_score2)

ssgsea_score2$sample <- gsub("\\.", "_", ssgsea_score2$sample)
ssgsea_score2$sample <- gsub("-", "_", ssgsea_score2$sample)

surv$sample <- rownames(surv)
surv$sample <- gsub("-", "_", surv$sample)


#整合表达谱数据和生存信息
train_v <- merge(ssgsea_score2, surv, by = "sample")

#删除没有生存信息的样本
table(is.na(train_v$OS))
table(train_v$OS, useNA = 'ifany')
#train_v <- subset(train_v, !is.na(train_v$OS))
table(is.na(train_v$OS.time))
#train_v <- subset(train_v, !is.na(train_v$OS.time))

#删除生存时间为0的样本
which(train_v$OS.time==0)
#train_v <- train_v[-262,]

#删除NA——没有检测到表达的基因
c <- which(grepl("NA",colnames(train_v)))
c
colnames(train_v)[c]

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

# save(diff_expr_clinical, file="diff_expr_clinical.rda")

covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))
univ_models <- lapply(univ_formulas,
                      function(x) coxph(x, data = diff_expr_clinical))

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
# write.csv(res, file = "03_uni_cox_OS_res.csv", quote = F)

# 选择单因素P值
res_results<- na.omit(res[which(as.numeric(res$p.value) < 0.05),])
rownames(res_results)
# write.csv(res_results, file = "04_uni_cox_OS_sig0.1.csv", quote = F)
dim(res_results)
# [1] 10 2

# PH检验
library(survminer)
library(gridExtra)
# RPUSD1 <- cox.zph(univ_models$RPUSD1)
# RPUSD1 <- ggcoxzph(RPUSD1)
# print(RPUSD1)
# ggsave("ggcoxzph_BLVRA.pdf",arrangeGrob(grobs = BLVRA))
# ggsave("ggcoxzph_BLVRA.png",arrangeGrob(grobs = BLVRA))
# + title = "BLVRA"
# ----------------------------
# PH检验
#rownames(res_results)
# which(colnames(df_merge)=="KAT2B") # 2
# which(colnames(df_merge)=="RPUSD1") # 13
# which(colnames(df_merge)=="RPLP0P6") # 20

x <- df_merge[,c(ncol(df_merge)-1,ncol(df_merge))]
uniSigExp = diff_expr_clinical
uniSigExp = uniSigExp[,colnames(uniSigExp) %in% rownames(res_results)]
uniSigExp <- merge(uniSigExp,x, by = "row.names")
rownames(uniSigExp) <- uniSigExp$Row.names
uniSigExp <- uniSigExp[,-1]
dat <- uniSigExp
outPH=data.frame()
for(i in colnames(dat[,1:(ncol(dat)-2)])){
  cox <- coxph(Surv(OS.time, OS) ~ dat[,i], data = dat)
  test.ph <- cox.zph(cox)
  #coxP=test.ph$coefficients[,"Pr(>|z|)"]
  outPH=rbind(outPH,
              cbind(id=i,
                    p=test.ph$table[1,"p"])
  )
}

sigPH=outPH[as.numeric(as.vector(outPH$p))>0.05,]
sigPH$id
write.csv(sigPH,file="03_PH.Sig.csv",quote=F)
sigPH <- sigPH[order(sigPH$p,decreasing = T),]
write.csv(sigPH, "03_sigPH.csv")

# 单因素cox森林图----
res_results<- res_results[order(res_results$p.value),]
# load("../01_DEG/RRGs.rda")
# intersect(rrgs, rownames(res_results))

res_results_plot <- res_results
res_results_plot <- res_results_plot[order(res_results_plot$p.value),]

unix_res <- tidyr::separate(res_results_plot, "HR (95% CI for HR)", into = c("HR", "HR.95L", "HR.95H"), sep = " ") %>% 
  tidyr::separate("HR.95L", into = c("HR.95L", "HR.95H"), sep = "\\-") %>% 
  tibble::rownames_to_column(var="GeneName")
unix_res$HR.95L <- gsub("\\(", "", unix_res$HR.95L)
unix_res$HR.95H <- gsub("\\)", "", unix_res$HR.95H)
unix_res[, 2:ncol(unix_res)] <- as.data.frame(apply(unix_res[, 2:ncol(unix_res)], 2, as.numeric))
unix_res <- unix_res[order(unix_res$p.value),]

hz <- paste(round(unix_res$HR,4),
            "(",round(unix_res$HR.95L,3),
            "-",round(unix_res$HR.95H,3),")",sep = "")
tabletext <- cbind(c(NA,"GeneName", unix_res$GeneName),
                   c(NA,"P value", ifelse(unix_res$p.value<0.001,
                                          "< 0.001",
                                          round(unix_res$p.value,4))),
                   c(NA,"Hazard Ratio(90% CI)",hz))

nrow(tabletext) + 1
tabletext[2,1] <- "Cell"

library(forestplot)
pdf(file = "05_univariate_cox_forest.pdf", height = 16, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE, rep(FALSE, 90)),
           col=fpColors(box="red", lines="royalblue", zero = "black"),
           mean=c(NA,NA,unix_res$HR),
           lower=c(NA,NA,unix_res$HR.95L), #95%置信区间下限
           upper=c(NA,NA,unix_res$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.2,"npc"), #图在表中的宽度比例
           cex=1.2, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontface = "bold", fontfamily="Times"),
                          summary = gpar(cex=1.2, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.2, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily="Times")),
           xlab="Hazard Ratio",
           cex.lab=1,
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()


diff_expr_clinical1 <- diff_expr_clinical
colnames(diff_expr_clinical1)
#年龄age
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "T_cell",  #需要计算的数据列名
                                   minprop  = 0.3
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res

diff_expr_clinical1$risk1 <- ifelse(diff_expr_clinical1[,"T_cell"] > res["T_cell",1],'High','Low')
table(diff_expr_clinical1$risk1)

library(survminer)
# K-M曲线 ----------
kmfit <- survfit(Surv(OS.time, OS) ~ risk1, data = diff_expr_clinical1)
surv_pvalue(kmfit)

pdf('04_KM_T_cell.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High","Low" ),
                     legend.title="Risk",
                     title="T_cell",
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

rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./04_GSE14520_SurvivalVerify")){
  dir.create("./04_GSE14520_SurvivalVerify")
}
setwd("./04_GSE14520_SurvivalVerify")

getwd()


data=read.table("../00.data/GEO_HCC/GSE14520/GSE14520_GPL3921.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
#intersect(colnames(data2),survival_data_val$sample)
survival_data_val=read.table("../00.data/GEO_HCC/GSE14520/survival.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
survival_data_val=survival_data_val[survival_data_val$sample %in% colnames(data),]
sample=intersect(colnames(data),survival_data_val$sample)
colnames(survival_data_val)=c("sample","OS","OS.time")
survival_data_val$OS.time <- survival_data_val$OS.time*30

str(survival_data_val)

#读取预后基因
lasso_geneids <- readRDS("../03_Survival/08_lasso_geneids.rds")


data1 <- data[,sample]
grep("FOXO1",rownames(data1))
rownames(data1)[3042] <- "FOXO1"

data1 <- data1[rownames(data1) %in% lasso_geneids,]

setdiff(lasso_geneids, rownames(data1))

verify_dat <- data1
verify_dat <- as.data.frame(t(verify_dat))

verify_dat <- rownames_to_column(as.data.frame(verify_dat),var = "sample")


surv <- survival_data_val
rownames(surv) <- surv$sample
identical(surv$sample,verify_dat$sample)
surv <- surv[verify_dat$sample,]
rt <- merge(surv,verify_dat,by = "sample") %>% na.omit()

rt$OS <- as.numeric(rt$OS)
rt$OS.time <- rt$OS.time/30

range(rt$OS.time)


save(rt, file = "01_verify_expr_survival.rda")
write.csv(rt, "01_verify_expr_survival.csv")

rt1 <- rt[,-c(1)]
rownames(rt1) <- rt$sample

rt2 <- rt1[,-c(1:2)]
rt2 <- as.data.frame(rt2)

# 计算风险评分------------------------------------------------------------------
# 行为样本，列为 OS.time和OS
risk <- rt1[,c(1,2)]
range(rt2)
# rt2 <- log2(rt2)
# risk$OS <- as.numeric(risk$OS)
# risk$OS.time <- as.numeric(risk$OS.time)
load("../03_Survival/10_Lasso_Coefficients.rda")
lasso <- df.coef
str(lasso)
lasso <- as.data.frame(lasso)
lasso$coef <- as.numeric(lasso$coefficient)
str(lasso$coef)
rt2 <- rt2[,lasso$gene]
cnt<-1
while (cnt < nrow(risk)+1) {
  risk$riskScore[cnt]<-sum(lasso$coef*subset(rt2,select=lasso$gene)[cnt,])
  cnt = cnt + 1
}
str(risk)
table(is.na(risk$riskScore))

#去除小于90天的样本
# risk <- risk[risk$OS.time>90,]

# 分组 ----------------------------------------------------------------------
# KM_group = "median"
KM_group = "cutpoint"
cut.num = 0.3
# if(KM_group == "median"){
#   
#   # 中位数分组 -------------------------------------------------------------------
#   risk$risk <- ifelse(risk$riskscore > median(risk$riskscore),"High risk","Low risk")
#   str(risk)
#   # write.table(risk,"04.risk.txt",sep='\t',quote=F)
#   saveRDS(risk,"04_risk_verify.rds")
# }else{
# 最佳截断值分组 -----------------------------------------------------------------
library(survival)
library(survminer)
res.cut <- survminer::surv_cutpoint(risk, #数据集
                                    minprop = cut.num,
                                    time = "OS.time", #生存时间
                                    event = "OS", #生存状态
                                    variables = "riskScore"  #需要计算的数据列名
)
res<-data.frame(summary(res.cut)) #查看数据最佳截断点及统计量
res # 0.860394
risk$risk <- ifelse(risk[,"riskScore"] > res["riskScore",1],'High','Low')


ROC <- risk
#生存时间使用天数，如果是月份需要乘以30,年份乘以365
range(ROC$OS.time)

cutoff_1 <- 1*12
cutoff_2 <- 2*12
cutoff_3 <- 3*12
cutoff_4 <- 4*12
cutoff_5 <- 5*12
cutoff_6 <- 6*12
cutoff_7 <- 7*12
cutoff_9 <- 9*12

year_1= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_1,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_2= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_2,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_3= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_3,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_4= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_4,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_5= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_5,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_6= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_6,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_7= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_7,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_9= survivalROC(Stime=ROC$OS.time,##生存时间
                    status=ROC$OS,## 终止事件
                    marker = ROC$riskScore, ## marker value
                    predict.time = cutoff_9,## 预测时间截点
                    method = 'KM')##span,NNE法的namda

year_1[["AUC"]] # 0.6465583
year_2[["AUC"]] # 0.6617107
year_3[["AUC"]] # 0.6652345
year_4[["AUC"]] # 0.6683036
year_5[["AUC"]] # 0.6626439
year_6[["AUC"]] # 
year_7[["AUC"]] # 
year_9[["AUC"]] # 

years <- c(1,2,3)

if(T){
  pdf(file = paste0("03_ModelROC_verify.pdf"),width = 7,height = 7)
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("03_ModelROC_verify.png"),width= 7, height= 7, units="in", res=300)
  dev.control("enable")
  par(mar = c(5,5,5,2))
  plot(year_1$FP, year_1$TP,
       type="l",col="red",xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Fraction",
       ylab="True Positive Fraction",
       main=paste0("Verify, Method = KM\n Year = ",years[1],",",years[2],",",years[3]),
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5
  )
  abline(0,1,col="gray",lty=2)
  lines(year_1$FP, year_1$TP, type="l",col="#BC3C29FF",xlim=c(0,1), ylim=c(0,1))
  lines(year_2$FP, year_2$TP, type="l",col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
  lines(year_3$FP, year_3$TP, type="l",col="#E18727FF",xlim=c(0,1), ylim=c(0,1))
  #lines(year_4$FP, year_4$TP, type="l",col="#B23AEE",xlim=c(0,1), ylim=c(0,1))
  #lines(year_5$FP, year_5$TP, type="l",col="#4C8045",xlim=c(0,1), ylim=c(0,1))
  legend(0.5,0.2,c(paste("AUC of ",years[1],"-year =",round(year_1$AUC,3)),
                   paste("AUC of ",years[2],"-year =",round(year_2$AUC,3)),
                   paste("AUC of ",years[3],"-year =",round(year_3$AUC,3))#,
                   #paste("AUC of 5 year =",round(year_4$AUC,3)),
                   #paste("AUC of 5 year =",round(year_5$AUC,3))
  ),
  x.intersp=1, y.intersp=0.8,
  lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF",'#E18727FF'),#,"#B23AEE","#FF0000"
  bty = "n",# bty框的类型
  seg.len=1,cex=1.2)#
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
}


##timeROC---------------------------------------
library(timeROC)
time_roc_res <- timeROC(
  T = ROC$OS.time,
  delta =ROC$OS,
  marker = ROC$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(years[1] * 12, years[2] * 12, years[3] * 12),
  ROC = TRUE,
  iid = TRUE
)

time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_2year = time_roc_res$TP[, 2],
  FP_2year = time_roc_res$FP[, 2],
  TP_3year = time_roc_res$TP[, 3],
  FP_3year = time_roc_res$FP[, 3]
  
)

library(ggplot2)
train_timeROC <- ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_2year, y = TP_2year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#E18727FF") +
  # geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "green") +
  # geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  ggtitle(paste0("Verify, Method = KM\n Year = ", years[1],",",years[2],",",years[3]))+
  annotate("text",family='Times',
           x = 0.75, y = 0.35, size = 6,
           label = paste0("AUC at ",years[1],"-year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",family='Times',
           x = 0.75, y = 0.25, size = 6,
           label = paste0("AUC at ",years[2],"-year = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",family='Times',
           x = 0.75, y = 0.15, size = 6,
           label = paste0("AUC at ",years[3],"-year = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  # annotate("text",family='Times',
  #          x = 0.75, y = 0.45, size = 6,
  #          label = paste0("AUC at 4 days = ", sprintf("%.2f", time_roc_res$AUC[[4]])), color = "green"
  # ) +
  # annotate("text",family='Times',
  #          x = 0.75, y = 0.55, size = 6,
  #          label = paste0("AUC at 5 days = ", sprintf("%.2f", time_roc_res$AUC[[5]])), color = "#E18"
  # ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",color = "black",size = 24,family='Times'),
    axis.text.x = element_text(color = "black",size = 19,family='Times'),
    axis.text.y = element_text(color = "black",size = 19,family='Times'),
    axis.title.x = element_text(face = "bold",color = "black",size = 22,family='Times'),
    axis.title.y = element_text(face = "bold",color = "black",size = 22,family='Times'),
  )
train_timeROC
ggsave('02_ModeltimeROC.pdf',train_timeROC,w=7,h=7)
ggsave('02_ModeltimeROC.png',train_timeROC,w=7,h=7)


save(risk, file = "verify_riskScore.rda")
write.csv(risk, "verify_riskScore.csv", quote = F)

# KM生存分析---------------------------------------------------------------------
library(survival)
library(survminer)
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = risk)
p1 <- surv_pvalue(kmfit, method = 'GB')
pvalue <- p1$pval.txt
pvalue

png('04_verify_KM.png',w=7,h=7,family='Times',units = 'in',res = 600)
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

pdf('04_verify_KM.pdf',w=7,h=7,onefile=F,family='Times')
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

# riskScore dis -----
rt <- risk
table(rt$risk)
lowrisk <- rt[which(rt$risk == "Low"),]
table(lowrisk$OS)

highrisk <- rt[-which(rt$risk == "Low"),]
table(highrisk$OS)

rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]

lowLength <- table(risk$risk)[2]
highLength <- table(risk$risk)[1]
color=as.vector(rt$OS)
# c("#FC8D62", "#66C2A5")
# "#0073C2FF","#A73030FF"
color[color==0]="#0072B5FF"
color[color==1]="#D62728"

pdf(file = paste0("05_verify_survStat.pdf"),width = 8,height = 6)
a <- dev.cur()
png(file = paste0("05_verify_survStat.png"),width= 8, height= 6, units="in", res=600)
dev.control("enable")
par(mar = c(5,5,5,2))
plot(rt$OS.time,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="OS time (Months)",
     col=color,
     ylim=c(0,max(rt$OS.time)*1.5),
     cex.axis = 1.2,
     cex.lab = 1.4,
)
abline(v=lowLength,lty=2)
legend(0,max(rt$OS.time)*1.4,c(paste("Dead"),
                               paste("Alive")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#D62728", "#0072B5FF"),
       bty = "n",
       seg.len=1,cex=1.5)#
title(main = "Verify",cex.main = 1.3)
dev.copy(which = a)
dev.off()
dev.off()
# 风险曲线图---------------------------------------------------------------------
line=rt[,"riskScore"]
#line[line>10]=10

pdf(file = paste0("06_riskScore.pdf"),width = 8,height = 6)
a <- dev.cur()
png(file = paste0("06_riskScore.png"),width= 8, height= 6, units="in", res=300)
dev.control("enable")
par(mar = c(10,5,5,2))

plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#6495ed",lowLength),
           rep("#dc143c",highLength)),
     cex.axis = 1.2,
     cex.lab = 1.4
)
legend("bottomright",max(risk$riskScore),c(paste("High"),
                                           paste("Low")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#dc143c","#6495ed"),
       bty = "n",
       seg.len=1,cex=1.5)
abline(h= res$cutpoint,v=lowLength,lty=2)
title(main = "Verify")

dev.copy(which = a)
dev.off()
dev.off()



risk <- risk[order(risk$risk,decreasing = T),]

expr <- subset(data, select=rownames(risk))
risk$sample <- rownames(risk)
a <- data.frame(row.names=risk$sample,Risk=risk$risk)

table(a$Risk)
gene_order <- c("RECQL4","CDKN2B","HMGB2","IGFBP3","SERPINE1",
                "IGF1","PPARGC1A","FOXO1")
select.expr <- expr[gene_order,]
ann.color <- list(Risk=c(High = "#dc143c",Low ="#6495ed"))
matrix<- t(scale(t(select.expr)))#归一化
matrix[matrix< -1] <- -1
matrix[matrix>1] <- 1

pheatmap::pheatmap(matrix,main= "Verify",
                   fontfamily= "Times",
                   fontsize=10,
                   annotation_colors=ann.color,
                   annotation_col=a ,
                   annotation_names_col = F,
                   annotation_names_row = F,
                   color = colorRampPalette(c("#3182BD", "#FFFFFF","#E6550D"))(500),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_colnames = F,
                   width=6,height=4,
                   file='07_heatmap.png')
pheatmap::pheatmap(matrix,main= "Verify" ,
                   fontfamily= "Times",
                   fontsize=10,
                   annotation_colors = ann.color,
                   annotation_col = a,
                   annotation_names_col = F,
                   annotation_names_row = F,
                   color = colorRampPalette(c("#3182BD", "#FFFFFF","#E6550D"))(500),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_colnames = F,
                   width=6,height=4,
                   file='07_heatmap.pdf')


save.image(file = "LIHC_GSE14520_Survival.RData")

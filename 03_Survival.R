rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./03_Survival")){
  dir.create("./03_Survival")
}
setwd("./03_Survival")

getwd()

df=read.table("../00.data/TCGA/TCGA-LIHC.survival.tsv",sep="\t",header=T,stringsAsFactors = F,quote="")
df=df[grep("-01",df$sample),]
df=df[,c(4,3,2)]
colnames(df)[1]="sample"
clinical=df[,c("sample","OS","OS.time")]
rownames(clinical)=clinical$sample
clinical=clinical[,-1]
##--mRNA--##
gene=read.table("../02.venn/sig_age_padj_freq3.xls",sep="\t",header=T,stringsAsFactors = F,quote="")[,1]
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
data=data[gene,]
colnames(data)=substring(colnames(data),1,12)
data=t(data)
sameSample=intersect(row.names(data),row.names(clinical))##365
data=data[sameSample,]
clinical=clinical[sameSample,]
out=cbind(clinical,data)
out=cbind(id=row.names(out),out)
colnames(out)[1] <- "sample"
write.table(out,file="OS_gene.txt",sep="\t",row.names=F,quote=F)

df_merge <- out

# unicox ------------------------------------------------------------------
library(survival)
diff_expr_clinical <- df_merge[,-1]
colnames_sum <- colnames(diff_expr_clinical)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(diff_expr_clinical) <- colnames_sum

str(diff_expr_clinical$OS.time)
table(diff_expr_clinical$OS)

save(diff_expr_clinical, file="diff_expr_clinical.rda")

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
write.table(res, file = "03_uni_cox_OS_res.txt", sep = "\t", quote = F, row.names = T)

# 选择单因素P值
res_results<- na.omit(res[which(as.numeric(res$p.value) < 0.05),])
rownames(res_results)
write.table(res_results, file = "04_uni_cox_OS_sig0.1.txt", sep = "\t", quote = F, row.names = T)
dim(res_results)


# PH检验
library(survminer)
library(gridExtra)
# ----------------------------
# PH检验
x <- diff_expr_clinical[,c(1,2)]
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
write.table(sigPH,file="05_PH.Sig.txt",sep="\t",row.names=F,quote=F)
sigPH <- sigPH[order(sigPH$p,decreasing = T),]
write.csv(sigPH, "05_sigPH.csv")
# -----------------------------------------

# 单因素cox森林图----
res_results<- res_results[order(res_results$p.value),]

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
                   c(NA,"Hazard Ratio(95% CI)",hz))
# tabletext[3,3] <- "0.000188 (7.39e-8-0.476)"
nrow(tabletext) + 1

library(forestplot)
pdf(file = "06_univariate_cox_forest.pdf", height = 16, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE, rep(FALSE, 54)),
           col=fpColors(box="red", lines="royalblue", zero = "black"),
           mean=c(NA,NA,unix_res$HR),
           lower=c(NA,NA,unix_res$HR.95L), #95%置信区间下限
           upper=c(NA,NA,unix_res$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 0.5, 1, 1.5, 2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
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
           cex.lab=2,
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()

png(filename = "06_univariate_cox_forest.png", height = 16, width = 10, units = "in", res = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE, rep(FALSE, 54)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,unix_res$HR),
           lower=c(NA,NA,unix_res$HR.95L), #95%置信区间下限
           upper=c(NA,NA,unix_res$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=1.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 0.5, 1, 1.5, 2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
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
           cex.lab=2,
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()


# 1.Lasso模型   ----------------
library(lance)
library(readxl)
library(readr)
library(tidyverse)


## 01.获取数据集 -----------------------------------------------------------  
#获取与生存显著相关的基因 单因素cox回归分析pvalue < 0.05，PH检验>0.05
cox_gene <- sigPH$id

train_data <- diff_expr_clinical
train_data <- train_data[,colnames(train_data) %in% cox_gene]

#合并生存数据
train_data$OS <- diff_expr_clinical$OS
train_data$OS.time <- diff_expr_clinical$OS.time

### 03.LASSO
library(survival)
library(survminer)
library(glmnet)


str(train_data$OS)
x_all <- subset(train_data, select = -c(OS, OS.time))
y_all <- subset(train_data, select = c(OS, OS.time))


##  04.拟合模型 ----------------------------------------------------------------
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox")

## 05.交叉验证 -----------------------------------------------------------------
set.seed(100)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=10,
                  family = "cox")#, alpha = 0.7
save(cvfit, file="07_lasso_result.rda")

plot(cvfit)

##画图
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1]  0.01866839

# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
length(coef.min)
coef.min


# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
# "RECQL4"   "CDKN2B"   "HMGB2"    "IGFBP3"   "IGF1"     "SERPINE1" "PPARGC1A" "FOXO1"


saveRDS(lasso_geneids, file= "08_lasso_geneids.rds")

lasso_geneids <- as.data.frame(lasso_geneids)
colnames(lasso_geneids) <- "symbol"

x <- coef(fit)
tmp <- as.data.frame(as.matrix(x)) 
tmp$coef <- row.names(tmp) 
tmp <- reshape::melt(tmp, id = "coef") 
tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
tmp$coef <- gsub('_','-',tmp$coef) 
tmp$lambda <- fit$lambda[tmp$variable+1] 

tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1]
# compute L1 norm
library(dplyr)
tmp <- tmp %>% filter(coef %in% lasso_geneids$symbol)
table(tmp$coef)

#图片美化
head(tmp)
pdf("09_lasso_model_name.pdf",height = 5, width = 8,family='Times')
ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Log(Lambda)") + 
  ylab('Coefficients')+ 
  theme_bw(base_rect_size = 2)+ 
  #scale_color_manual(values= c('#ff9898','#dedb8e','#99e0ab','#D94F04','#007172',
  #                             '#025259','#c49d93','#aec6e8','#F2C6C2','#86A69D'))+ 
  scale_x_continuous(expand = c(0.01,0.01))+ 
  scale_y_continuous(expand = c(0.01,0.01))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=15,color='black'), 
        axis.text = element_text(size=12,color='black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size=12,color='black'), 
        legend.position = 'right')+ 
  annotate('text',x = -4,y=-0.25,
           label=paste0('Optimal Log(Lambda) = ', log(cvfit$lambda.min)),
           color='black')+ 
  guides(col=guide_legend(ncol = 1))
dev.off()


png("09_lasso_model_name.png", height = 5, width = 8, units = 'in',res = 600,family='Times')
ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Log(Lambda)") + 
  ylab('Coefficients')+ 
  theme_bw(base_rect_size = 2)+ 
  #scale_color_manual(values= c('#ff9898','#dedb8e','#99e0ab','#D94F04','#007172',
  #                             '#025259','#c49d93','#aec6e8','#F2C6C2','#86A69D'))+ 
  scale_x_continuous(expand = c(0.01,0.01))+ 
  scale_y_continuous(expand = c(0.01,0.01))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=15,color='black'), 
        axis.text = element_text(size=12,color='black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size=12,color='black'), 
        legend.position = 'right')+ 
  annotate('text',x = -4,y=-0.25,
           label=paste0('Optimal Log(Lambda) = ', log(cvfit$lambda.min)),
           color='black')+ 
  guides(col=guide_legend(ncol = 5))
dev.off()


pdf("09_1_lasso_model_name.pdf", height = 6, width = 8)
plot(cvfit$glmnet.fit, xvar = 'lambda')
dev.off()

png("09_1_lasso_model_name.png", height = 6, width = 8, units = 'in',res = 600,family='Times')
plot(cvfit$glmnet.fit, xvar = 'lambda')
dev.off()

pdf("09_lasso_coef.pdf", width = 8, height = 6)
plot(cvfit,las=1)
text(x = (log(cvfit$lambda.min)-1), y = 100, labels = paste("log(lambda.min) =", log(cvfit$lambda.min)), pos = 4)
dev.off()

png("09_lasso_coef.png", width = 8, height = 6, units = 'in',res = 600)
plot(cvfit,las=1)
text(x = (log(cvfit$lambda.min)-1), y = 100, labels = paste("log(lambda.min) =", log(cvfit$lambda.min)), pos = 4)
dev.off()

#交叉验证图
lambda_vals <- cvfit$lambda
cv_means <- cvfit$cvm
cv_sds <- cvfit$cvsd
cv_up <- cv_means + cv_sds
cv_lo <- cv_means - cv_sds
nonzero <- cvfit$nzero

# Data frame for plotting
cv_data <- data.frame(log_lambda = log(lambda_vals), 
                      mean = cv_means, 
                      sd = cv_sds,
                      upper = cv_up,
                      lower = cv_lo,
                      nonzero = nonzero)

# Convert nonzero to factor for coloring
cv_data$nonzero <- as.factor(cv_data$nonzero)
levels(cv_data$nonzero) <- paste0(levels(cv_data$nonzero), " vars")

# Plotting
library(ggsci)
lamda_p <- ggplot(cv_data, aes(x = log_lambda, y = mean, color = nonzero)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, size = 1, alpha = 0.6) +
  geom_point(size = 2) +
  geom_vline(xintercept = log(cvfit$lambda.min), color = "red", linetype = "dashed", size = 0.8) +
  #geom_vline(xintercept = log(res.lasso$lambda.1se), color = "green", linetype = "dashed", size = 0.8) +
  xlab("Log Lambda") +
  ylab("Binomial Deviance") +
  ggtitle("Cross-Validation Plot") +
  theme_bw(base_rect_size = 1.5) +
  scale_color_manual(values = c(pal_npg()(10), pal_d3()(10), pal_lancet()(10), pal_aaas()(10))) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size = 15, color = 'black'), 
        axis.text = element_text(size = 12, color = 'black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, color = 'black'), 
        legend.position = 'right') + 
  # annotate('text', x = log(lasso_fit$lambda.min) + 2, y = max(cv_data$mean) - 0.1, 
  #          label = paste0("Optimal Lambda = ", round(lasso_fit$lambda.min, 10)), color = 'black') +
  guides(col = guide_legend(ncol = 1))+
  annotate('text',x = -4,y=12,
           label=paste0('Optimal Log(Lambda) = ', log(cvfit$lambda.min)),
           color='black')
lamda_p
# Save the plot
# ggsave("cross_validation_plot.png", width = 6, height = 5, dpi = 600, units = "in", bg = "white")
ggsave("09_1_lasso_coef.png", lamda_p, width = 8, height = 6, dpi = 600, units = "in", bg = "white")
ggsave("09_1_lasso_coef.pdf", lamda_p, width = 10, height = 6)

#--------------------------- 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = cvfit$lambda.min)  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min

df.coef = cbind(gene = rownames(coef.min), coefficient = coef.min[,1]) %>% as.data.frame()
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame
write.table(df.coef, "10_Lasso_Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)  #6个
save(df.coef, file = "10_Lasso_Coefficients.rda")

# -----------------------训练集验证-----------
diff_expr_clinical1 <- train_data
# riskScore <- predict(step_Cox,type="lp",newdata=diff_expr_clinical)
riskScore <- predict(cvfit, newx = as.matrix(x_all), s=cvfit$lambda.min)
riskScore <- data.frame(riskScore)

identical(rownames(riskScore), rownames(diff_expr_clinical1))
diff_expr_clinical1$riskScore <- riskScore[,1]

# 中位数
# diff_expr_clinical1$risk <- ifelse(diff_expr_clinical1$riskScore>median(diff_expr_clinical1$riskScore), "High", "Low")
# 最佳截断值
OS.cut <- survminer::surv_cutpoint(diff_expr_clinical1, #数据集
                                   minprop = 0.3,
                                   time = "OS.time", #生存时间
                                   event = "OS", #生存状态
                                   variables = "riskScore"  #需要计算的数据列名
)
res<-data.frame(summary(OS.cut)) #查看数据最佳截断点及统计量
res # 1.139719
diff_expr_clinical1$risk <- ifelse(diff_expr_clinical1[,"riskScore"] > res["riskScore",1],'High','Low')


write.table(diff_expr_clinical1, file = "11_all_train_riskScore.txt", sep = "\t", row.names = T, col.names = T, quote = F)
save(diff_expr_clinical1, file = "11_all_train_riskScore.rda")
#diff_expr_clinical <- read.table("10_all_train_riskScore.txt")

diff_expr_clinical1$OS.time <- diff_expr_clinical1$OS.time/30

kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = diff_expr_clinical1)
surv_pvalue(kmfit, method = 'GB')
surv_pvalue(kmfit, method = 'log-rank')



png('12_train_KM.png',w=7,h=7,family='Times',units = 'in',res = 600)
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

pdf('12_train_KM.pdf',w=7,h=7,onefile=F,family='Times')
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

# ROC----
library(survivalROC)
rt = subset(diff_expr_clinical1, select = c(OS, OS.time, riskScore))

ROC <- rt
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

year_1[["AUC"]] # 0.7376049
year_2[["AUC"]] # 0.6870299
year_3[["AUC"]] # 0.7027867
year_4[["AUC"]] # 0.7095097
year_5[["AUC"]] # 0.6433623
year_6[["AUC"]] # 0.6065186
year_7[["AUC"]] # 0.7122132
year_9[["AUC"]] # 0.766889
summary(year_1$TP)

if(T){
  pdf(file = paste0("13_ModelROC_train.pdf"),width = 8,height = 8)
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("13_ModelROC_train.png"),width= 8, height= 8, units="in", res=600)
  dev.control("enable")
  par(mar = c(5,5,5,2))
  plot(year_1$FP, year_1$TP,
       type="l",col="red",xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Fraction",
       ylab="True Positive Fraction",
       main="LIHC-train, Method = KM\n Year = 1,2,3",
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5
  )
  abline(0,1,col="gray",lty=2)
  lines(year_1$FP, year_1$TP, type="l",col="#BC3C29FF",xlim=c(0,1), ylim=c(0,1))
  lines(year_2$FP, year_2$TP, type="l",col="#0072B5FF",xlim=c(0,1), ylim=c(0,1))
  lines(year_3$FP, year_3$TP, type="l",col="#4C8045",xlim=c(0,1), ylim=c(0,1))
  # lines(year_4$FP, year_4$TP, type="l",col="#B23AEE",xlim=c(0,1), ylim=c(0,1))
  # lines(year_5$FP, year_5$TP, type="l",col="#FF0000",xlim=c(0,1), ylim=c(0,1))
  legend(0.6,0.2,c(paste("AUC of 1-year =",round(year_1$AUC,3)),
                   paste("AUC of 2-year =",round(year_2$AUC,3)),
                   paste("AUC of 3-year =",round(year_3$AUC,3))#,
                   # paste("AUC of 4 year =",round(year_4$AUC,3)),
                   # paste("AUC of 5 year =",round(year_5$AUC,3))
  ),
  x.intersp=1, y.intersp=0.8,
  lty= 1 ,lwd= 2,col=c("#BC3C29FF","#0072B5FF",'#E18727FF'),#,"#B23AEE","#FF0000"
  bty = "n",# bty框的类型
  seg.len=1,cex=1.2)#
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
}


library(timeROC)
time_roc_res <- timeROC(
  T = ROC$OS.time,
  delta =ROC$OS,
  marker = ROC$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(1 * 12, 2 * 12, 3 * 12),
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
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_2year, y = TP_2year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#E18727FF") +
  # geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "green") +
  # geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  ggtitle(paste0("LIHC ",', Method= KM \n Year = 1,2,3'))+
  annotate("text",family='Times',
           x = 0.75, y = 0.35, size = 6,
           label = paste0("AUC at 1-year = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",family='Times',
           x = 0.75, y = 0.25, size = 6,
           label = paste0("AUC at 2-year = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",family='Times',
           x = 0.75, y = 0.15, size = 6,
           label = paste0("AUC at 3-year = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
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

ggsave('13.Model_timeROC.pdf',w=8,h=8)
ggsave('13.Model_timeROC.png',w=8,h=8)

# riskScore dis -----
rt <- subset(diff_expr_clinical1, select = c(OS, OS.time, riskScore, risk))
table(rt$risk)
lowrisk <- rt[which(rt$risk == "Low"),]
table(lowrisk$OS)

highrisk <- rt[-which(rt$risk == "Low"),]
table(highrisk$OS)


rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]

lowLength <- table(rt$risk)[2]
highLength <- table(rt$risk)[1]
color=as.vector(rt$OS)
# c("#FC8D62", "#66C2A5")
# "#0073C2FF","#A73030FF"
color[color==0]="#0073C2FF"
color[color==1]="#A73030FF"

pdf(file = paste0("14_train_survStat.pdf"),width = 8,height = 5)
a <- dev.cur()
png(file = paste0("14_train_survStat.png"),width= 8, height= 5, units="in", res=300)
dev.control("enable")
par(mar = c(5,5,5,2))
plot(rt$OS.time,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="OS time (Months)",
     col=color,
     ylim=c(0,max(rt$OS.time)*2),
     cex.axis = 1.2,
     cex.lab = 1.4,
)
abline(v=lowLength,lty=2)
legend(0,max(rt$OS.time)*1.8,c(paste("Dead"),
                               paste("Alive")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#A73030FF", "#0073C2FF"),
       bty = "n",
       seg.len=1,cex=1.5)
title(main = "LIHC-Train",cex.main = 1.3)
dev.copy(which = a)
dev.off()
dev.off()
# 风险曲线图---------------------------------------------------------------------
# rt<-risk
# lowrisk <- rt[which(rt$risk == "Low risk"),]
# table(lowrisk$fustat)
# highrisk <- rt[-which(rt$risk == "Low risk"),]
# table(highrisk$fustat)
# 
# rt=rt[order(rt$riskscore),]
# riskClass=rt[,"risk"]
# 
# highLength=length(riskClass[riskClass=="High risk"])
# lowLength=length(riskClass[riskClass=="Low risk"])
line=rt[,"riskScore"]
#line[line>10]=10

pdf(file = paste0("15_riskScore.pdf"),width = 8,height = 5)
a <- dev.cur()
png(file = paste0("15_riskScore.png"),width= 8, height= 5, units="in", res=300)
dev.control("enable")
par(mar = c(10,5,5,2))

plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#0073C2FF",lowLength),
           rep("#A73030FF",highLength)),
     cex.axis = 1.2,
     cex.lab = 1.4,
)
legend("bottomright",max(rt$riskScore),c(paste("High"),
                                         paste("Low")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#A73030FF","#0073C2FF"),
       bty = "n",
       seg.len=1,cex=1.5)#"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"
abline(h= res$cutpoint,v=lowLength,lty=2)
title(main = "LIHC-Train")

dev.copy(which = a)
dev.off()
dev.off()  

#symbol <- lasso_geneids$symbol
save(lasso_geneids, file="LASSO_GeneID.rda")


risk <- ROC
# risk$sample <- df_merge$sample
save(risk, file = "train_riskScore.rda")

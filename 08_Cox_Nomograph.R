rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./08_Cox_Nomograph")){
  dir.create("./08_Cox_Nomograph")
}
setwd("./08_Cox_Nomograph")

getwd()

#读取临床信息
# clinical <- read.table("../00.data/TCGA/TCGA-LIHC-clinical.txt", header = T, sep = "\t")
clinical <- read.delim2(gzfile("/data/nas2/database/TCGA/TCGA-LIHC.clinical.tsv.gz"),
                        header = T, check.names = F)
#读取风险得分
load("../03_Survival/11_all_train_riskScore.rda")


#提取样本
clinical1 <- clinical[clinical$submitter_id %in% rownames(diff_expr_clinical1),]

colnames(clinical1)


lihc_clinical <- data.frame(sample=clinical1$submitter_id,
                            age=clinical1$age_at_index,
                            gender=clinical1$gender,
                            weight=clinical1$initial_weight.samples,
                            race=clinical1$race,
                            Stage=clinical1$ajcc_pathologic_stage
                            # T_stage=clinical1$ajcc_pathologic_t,
                            # M_stage=clinical1$ajcc_pathologic_m,
                            # malignancy=clinical1$prior_malignancy
) %>% dplyr::distinct(sample, .keep_all = T)


#读取风险评分
risk <- diff_expr_clinical1[,c("OS","OS.time","riskScore")]
risk$sample <- rownames(risk)


lihc_clinical <- merge(lihc_clinical, risk, by = "sample")

save(lihc_clinical, file = "01_lihc_clinical_infor.rda")
write.csv(lihc_clinical, "01_lihc_clinical_infor.csv")

# 独立预后----
library(magrittr)
df_all <- lihc_clinical
str(df_all)

df_all[df_all == ""] <- NA
df_all[df_all == "NA"] <- NA


colnames(df_all)

table(df_all$gender, useNA="ifany")
df_all$gender <- as.factor(df_all$gender)

table(df_all$race, useNA="ifany")
df_all$race <- gsub("not reported", NA, df_all$race)
df_all$race <- factor(df_all$race,
                      levels = c("american indian or alaska native","asian",
                                 "black or african american","white"),
                      labels = c("Others","Others",
                                 "Others","White"
                      ))

# table(df_all$T_stage, useNA="ifany")
# df_all$T_stage <- gsub("TX", NA, df_all$T_stage)
# df_all$T_stage <- factor(df_all$T_stage,
#                          levels = c("T1","T2","T2a","T2b",
#                                     "T3","T3a","T3b","T4"),
#                          labels = c("T1","T2|3|4","T2|3|4","T2|3|4",
#                                     "T2|3|4","T2|3|4","T2|3|4","T2|3|4"
#                                    ))

# table(df_all$M_stage, useNA="ifany")
# df_all$M_stage <- gsub("TX", NA, df_all$M_stage)
# df_all$M_stage <- factor(df_all$M_stage,
#                          levels = c("M0","M1"),
#                          labels = c("M0","M1"
#                          ))

table(df_all$Stage, useNA="ifany")
# df_all$Stage <- gsub("not reported", NA, df_all$Stage)
df_all$Stage <- factor(df_all$Stage,
                       levels = c("Stage I","Stage II",
                                  "Stage III","Stage IIIA","Stage IIIB","Stage IIIC",
                                  "Stage IV","Stage IVA","Stage IVB"),
                       labels = c("Stage1","Stage2",
                                  "Stage3|4","Stage3|4","Stage3|4","Stage3|4",
                                  "Stage3|4","Stage3|4","Stage3|4"))


# table(df_all$malignancy, useNA="ifany")
# df_all$malignancy <- as.factor(df_all$malignancy)


df_all$age <- as.numeric(df_all$age)
range(df_all$age)
df_all$age <- ifelse(df_all$age>=70,">=70","<70")
df_all$age <- factor(df_all$age)

df_all$weight <- as.numeric(df_all$weight)
range(df_all$weight)

str(df_all)


library(survival)
str(df_all)
# rownames(df_all) <- df_all$sample
# df_all <- df_all[,-1]
save(df_all, file = "02_train_clinical_character_cox.rda")

coxvariates <- colnames(df_all)[!colnames(df_all) %in% c("sample", "OS", "OS.time")]

uncox.formula <- sapply(coxvariates, function(x) {as.formula(paste("Surv(OS.time, OS)~", x))})

res.ref = c(1,1,1,NA)
univ_models <- list()  # 新增：保存所有模型
table2 <- NULL
merge <- NULL
for (i in 1:length(coxvariates)){
  name2 <- coxvariates[i]
  # formula <- as.formula(paste0("Surv(time = OS.time, event = OS) ~ `", name2, "`"))
  # model <- coxph(formula, data = os_risk_clinical)
  
  res = coxph(as.formula(paste0("Surv(time = OS.time, event = OS) ~ `",name2,"`")), data = df_all)
  res2 <- res
  res <- res%>% summary
  
  univ_models[[name2]] <- res2
  if(length(table(df_all[[name2]]))==2|name2%in%c("lymph_node",'age',"weight",'riskScore')){
    res = c(res$conf.int[-2], res$coefficients[5])
  }
  if(length(table(df_all[[name2]]))>2 & !name2%in%c("lymph_node",'age',"weight",'riskScore')){
    res = cbind(res$conf.int[,-2], res$coefficients[,5])
  }
  ph <- cox.zph(coxph(as.formula(paste0("Surv(time = OS.time,event = OS) ~ `",name2,"`")),data=df_all))$table[1,,drop = F]
  res2 <- ph
  merge1 <- res
  p <- ph%>%as.data.frame()
  
  
  if(is.null(table2)){
    table2 <- res2
  }else{
    table2 <- rbind(table2,res2)
  }
  if(p$p >= 0.05){
    if(length(table(df_all[[name2]]))==2|name2%in%c("lymph_node",'age',"weight",'riskScore')){
      merge1 <-  matrix(merge1, nrow = 1)
      rownames(merge1) <- name2
    }
    if(length(table(df_all[[name2]]))>2 & !name2%in%c("lymph_node",'age',"weight",'riskScore')){
      merge1 <- rbind(res.ref,merge1)
      rownames(merge1) <- c(paste0(name2,'(Reference)'),names(table(df_all[[name2]]))[-1])
    }
    
    
    if(is.null(merge)){
      merge <- merge1
    }else{
      merge <- rbind(merge1,merge)
    }
  }
}
print(table2)
write.csv(table2,file = '03_univariate_cox_ph.csv')

library(ggplot2)
library(cowplot)

ph.plot <- lapply(univ_models, function(model){
  zph.test <- cox.zph(model)
  plot_list <- ggcoxzph(zph.test)
  plot_list[[1]]  # 提取主图
})

combined.plot <- plot_grid(plotlist = ph.plot, nrow = 3)
combined.plot
# 保存图像
ggsave(filename = "04_uncox_ph.pdf", combined.plot, width = 14.5, height = 10)
ggsave(filename = "04_uncox_ph.png", combined.plot, width = 14.5, height = 10, bg = "white", units = "in")


rownames(merge)
res <- merge%>%as.data.frame()
res$Indicators <- rownames(res)

colnames(res) = c("hr","low","up","pv","Indicator")
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
res2 <- subset(res2, select = -c(Indicator))

write.csv(res2, file = "05_univariate_cox_prog_forest.csv", quote = F, row.names = T)

#整理画图数据
library(tidyr)
#整理成HR(HR.95L-HR.95H)
hz <- paste(round(res2$HR,4), #保留3位小数
            "(",
            round(res2$HR.95L,4),
            "-",
            round(res2$HR.95H,4),
            ")",
            sep = "") 
hz

# hz[c(4)] <- "" #对照值为空

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.05,"< 0.05",round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
tabletext
nrow(tabletext) + 1

#画图
library(forestplot)
pdf(file = '06_uncoxforest.pdf',height = 10, width = 14, onefile = F,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(FALSE,rep(FALSE,1),rep(FALSE, 2),rep(F,1),rep(FALSE,12)),
           col=fpColors(box="#D62728", lines="#0072B5FF", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=1.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,2,4,6), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2, fontfamily="Times"),
                          summary = gpar(cex=1.4, fontfamily = "Times"),
                          ticks=gpar(cex=1.2, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.3, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.5, fontface = "bold", fontfamily="Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate Cox",
           clip = c(0,2)) # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '06_uncoxforest.png',height = 8, width = 12, units = 'in',res = 600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(FALSE,rep(FALSE,1),rep(FALSE, 2),rep(F,1),rep(FALSE,12)),
           col=fpColors(box="#D62728", lines="#0072B5FF", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=1.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,2,4,6), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1.2, fontfamily="Times"),
                          summary = gpar(cex=1.4, fontfamily = "Times"),
                          ticks=gpar(cex=1.2, fontface = "bold", fontfamily="Times"),
                          xlab=gpar(cex = 1.3, fontface = "bold", fontfamily="Times"),
                          title=gpar(cex = 1.5, fontface = "bold", fontfamily="Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate Cox",
           clip = c(0,2)) # 垂直于x轴的网格线，对应每个刻度
dev.off()


# 多因素----
df_all1 <- df_all
colnames(df_all1)
features <- c("age","weight","Stage","riskScore")
cox_data <- as.formula(paste0('Surv(OS.time, OS)~', paste(features, collapse = "+")))
cox_more <- coxph(cox_data, data = df_all1)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]

cox_formula <- as.formula(paste("Surv(OS.time, OS)~",
                                paste(rownames(cox_table),
                                      collapse = "+")))
cox_more_2 <- coxph(cox_formula, data = df_all1)
mul_cox_result <- summary(cox_more_2)$coefficients
mul_cox_cli <- rownames(mul_cox_result)
mul_cox_cli
library(survminer)
# library(ezcox)
# show_models(cox_more_2)
pdf(file = "07_multivariate_cox_prog_forest.pdf", family = "Times", height = 4, width = 8, onefile = F)
ggforest(model = cox_more_2,
         data = df_all1,
         main = "Hazard ratio of Candidate Feature",
         fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
dev.off()
png(file = "07_multivariate_cox_prog_forest.png", family = "Times", height = 4, width = 8, units = "in", res = 600)
ggforest(model = cox_more_2,
         data = df_all1,
         main = "Hazard ratio of Candidate Feature",
         fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
dev.off()

#重新绘制多因素森林图
res.mul = coxph(cox_formula, data = df_all1) %>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
res.mul
res.mul = rbind(res.ref,res.mul[1,],
                res.mul[2,],
                res.ref,res.mul[3:4,],
                res.mul[5,]
)

res.mul
res.mul$Indicators = c('age(age<70)','age(age>=70)',
                       'weight',
                       'Stage(Stage1)','Stage2','Stage3|4',
                       'riskScore'
)
colnames(res.mul) = c("hr","low","up","pv","Indicator")

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)

rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "08_multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,7),
            "(",round(multi_res$HR.95L,7),
            "-",round(multi_res$HR.95H,7),")",sep = "")
hz

tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.05,
                                      "< 0.05",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
tabletext
nrow(tabletext) + 1
library(forestplot)
pdf(file = "09_multivariate_cox_prog_forest.pdf", family = "Times", height = 7, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,rep(FALSE,10),TRUE,rep(FALSE,10),TRUE,rep(FALSE,10),TRUE,rep(FALSE,10)),
           col=fpColors(box="#D62728", lines="#0072B5FF", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=2,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2, 3, 4, 5), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily = "Times"),
                          xlab=gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate",
           clip = c(0,30)) # 垂直于x轴的网格线，对应每个刻度

dev.off()
png(filename = "09_multivariate_cox_prog_forest.png", family = "Times", height = 7, width = 12, units = "in", res = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,rep(FALSE,10),TRUE,rep(FALSE,12),TRUE,rep(FALSE,14)),
           col=fpColors(box="#D62728", lines="#0072B5FF", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=2,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2, 3, 4, 5), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily = "Times"),
                          xlab=gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate",
           clip = c(0,30)) # 垂直于x轴的网格线，对应每个刻度

dev.off()

save(cox_more_2, file = "10_multicox.rda")
save(cox_zph, file = "11_multicox_phtest.rda")
cox_zph

# 生成Schoenfeld残差图列表
ph_plots <- ggcoxzph(cox_zph)

# 合并所有残差图，2行排列
combined_plot <- plot_grid(plotlist = ph_plots, nrow = 2)+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 16))

# 保存图形文件
ggsave(filename = "11_multi_cox_ph.pdf", combined_plot, width = 10, height = 8)
ggsave(filename = "11_multi_cox_ph.png", combined_plot, width = 10, height = 8, bg = "white", units = "in")


#列线图模型
#选取单因素和多因素都通过的临床特征 age,path_N,riskScore
colnames(df_all)
df_all2 <- df_all


# Nomogram----
library(rms)
ddist <- datadist(df_all2)
# df_all2$OS.time <- df_all2$OS.time*365
options(datadist = "ddist")
cox_data2 <- as.formula(paste0('Surv(OS.time, OS)~', paste(c("age","weight","Stage","riskScore"), collapse = "+")))
res.cox <- psm(cox_data2, data = df_all2, dist = "lognormal")
surv <- Survival(res.cox)

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(1*12, x), function(x) surv(2*12, x), function(x) surv(3*12, x)), #记得修改时间，年，月，日
                    funlabel = c("1-year Survival Probability", "2-year Survival Probability", "3-year Survival Probability"),
                    lp = F,
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99)
)


library(regplot)

f <- cph(Surv(OS.time, OS)~ riskScore + Stage + weight + age,data = df_all2,x =T,y=T,surv = T)
f
surv <- Survival(f)
surv1 <- function(x) surv(1*12,x)
surv2 <- function(x) surv(2*12,x)
surv3 <- function(x) surv(3*12,x)

regplot(f,failtime = c(1*12,3*12,5*12), plots=c("violin","bars"), showP = T,
        observation=TRUE, title="CRCLM Survial Nomogram", clickable=TRUE, points=TRUE,droplines=TRUE)

regplot(f, failtime = c(1*12,2*12,3*12), plots=c("violin","bars"), observation=df_all[10,], 
        title="Survial Nomogram",prfail=T, showP = T, points=TRUE, droplines=TRUE)
my_plot <- recordPlot()  # 保存当前绘图
png("12_Nomograph.png",width=10,height=8, units = "in", res = 600, family = "Times")
replayPlot(my_plot)
dev.off()
pdf("12_Nomograph.pdf",width=10,height=8, family = "Times")
replayPlot(my_plot)
dev.off()

set.seed(123)

coxfit1 <- cph(Surv(OS.time, OS)~ age + weight + Stage + riskScore,data = df_all2, x=T,y=T,surv = T,
               time.inc = 1*12
)
# 参数u需要改成自己使用的年份或者月份或者天数
cal1 <- rms::calibrate(coxfit1, cmethod="KM", method="boot",u=1*12,m=50,B=500)


coxfit2 <- cph(Surv(OS.time, OS)~ age + weight + Stage + riskScore,data = df_all2, x=T,y=T,surv = T,
               time.inc = 2*12
)

cal2 <- rms::calibrate(coxfit2, cmethod="KM", method="boot",u=2*12,m=50,B=500)


coxfit3 <- cph(Surv(OS.time, OS)~ age + weight + Stage + riskScore,data = df_all2, x=T,y=T,surv = T,
               time.inc = 3*12
)

cal3 <- rms::calibrate(coxfit3, cmethod="KM", method="boot",u=3*12,m=50,B=500)

png(filename = "13_nomogram_predicted.png", family = "Times", height = 8, width = 8, units = "in", res = 600)
par(mar=c(5,4,2,3),cex=1.5,family="Times")
plot(cal1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal2,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
#加上图例
legend("bottomright", legend=c("1-year", "2-year", "3-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

pdf(file = "13_nomogram_predicted.pdf", family = "Times", height = 8, width = 8)
par(mar=c(5,4,2,3),cex=1.5,family="Times")
plot(cal1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal2,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
#加上图例
legend("bottomright", legend=c("1-year", "2-year", "3-year"),
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

#  ROC----
library(survivalROC)
library(tidyverse)

rt = subset(df_all2, select = c(OS, OS.time, riskScore))
cutoff_1 <- 1*12 #修改年份，月份，天数
cutoff_2 <- 2*12
cutoff_3 <- 3*12

year_1= survivalROC(Stime=rt$OS.time,##生存时间
                    status=rt$OS,## 终止事件
                    marker = predict(f, type = "lp"), ## marker value
                    predict.time = cutoff_1,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_2= survivalROC(Stime=rt$OS.time,##生存时间
                    status=rt$OS,## 终止事件
                    marker = predict(f, type = "lp"), ## marker value
                    predict.time = cutoff_2,## 预测时间截点
                    method = 'KM')##span,NNE法的namda
year_3= survivalROC(Stime=rt$OS.time,##生存时间
                    status=rt$OS,## 终止事件
                    marker = predict(f, type = "lp"), ## marker value
                    predict.time = cutoff_3,## 预测时间截点
                    method = 'KM')##span,NNE法的namda

if(T){
  pdf(file = paste0("14_ModelROC_train.pdf"),width = 7,height = 7)
  a <- dev.cur()   #记录pdf设备
  png(file = paste0("14_ModelROC_train.png"),width= 7, height= 7, units="in", res=600)
  dev.control("enable")
  par(mar = c(5,5,5,2))
  plot(year_1$FP, year_1$TP,
       type="l",col="red",xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Fraction",
       ylab="True Positive Fraction",
       main="Nomogram, Method = KM\n Year = 1,2,3",  #修改疾病和年份
       cex.lab = 1.5,
       cex.axis = 1.5,
       cex.main = 1.5
  )
  abline(0,1,col="gray",lty=2)
  lines(year_1$FP, year_1$TP, type="l",col="#EB4B17",xlim=c(0,1), ylim=c(0,1))
  lines(year_2$FP, year_2$TP, type="l",col="#2775AB",xlim=c(0,1), ylim=c(0,1))
  lines(year_3$FP, year_3$TP, type="l",col="#4C8045",xlim=c(0,1), ylim=c(0,1))
  # lines(year_4$FP, year_4$TP, type="l",col="#B23AEE",xlim=c(0,1), ylim=c(0,1))
  # lines(year_5$FP, year_5$TP, type="l",col="#FF0000",xlim=c(0,1), ylim=c(0,1))
  legend(0.5,0.2,c(paste("AUC of 1-year =",round(year_1$AUC,3)),
                   paste("AUC of 2-year =",round(year_2$AUC,3)),
                   paste("AUC of 3-year =",round(year_3$AUC,3))#,
                   # paste("AUC of 4 year =",round(year_4$AUC,3)),
                   # paste("AUC of 5 year =",round(year_5$AUC,3))
  ),
  x.intersp=1, y.intersp=0.8,
  lty= 1 ,lwd= 2,col=c("#EB4B17","#2775AB",'#4C8045'),#,"#B23AEE","#FF0000"
  bty = "n",# bty框的类型
  seg.len=1,cex=1.2)#
  dev.copy(which = a)  #复制来自png设备的图片到pdf
  dev.off()
  dev.off()
}


years <- c(1,2,3)

##timeROC---------------------------------------
library(timeROC)
colnames(df_all2)
res.mul1 <- coxph(cox_data2, data = df_all2)
identical(df_all2$sample, risk$sample)
df_all3 <- df_all2[match(risk$sample, df_all2$sample),]
identical(df_all3$sample, risk$sample)

risk$pred <- predict(res.mul1, newdata = df_all3, type = "lp")#计算列线图得分
cli_dat <- risk

# cli_dat$OS.time <- cli_dat$OS.time/365
colnames(cli_dat)[colnames(cli_dat)=="pred"] <- "nomogram"

ROC.nomo <- timeROC(T=cli_dat$OS.time,
                    delta=cli_dat$OS,
                    marker=cli_dat$nomogram,
                    cause=1,
                    weighting="marginal",
                    times=c(years[1] * 12, years[2] * 12, years[3] * 12),
                    ROC = TRUE,
                    iid=TRUE)

time_ROC_df <- data.frame(
  TP_1year = ROC.nomo$TP[, 1],
  FP_1year = ROC.nomo$FP[, 1],
  TP_2year = ROC.nomo$TP[, 2],
  FP_2year = ROC.nomo$FP[, 2],
  TP_3year = ROC.nomo$TP[, 3],
  FP_3year = ROC.nomo$FP[, 3]
  
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
  ggtitle(paste0("Nomogram, Method = KM\n Year = ", years[1],",",years[2],",",years[3]))+
  annotate("text",family='Times',
           x = 0.75, y = 0.35, size = 6,
           label = paste0("AUC at ",years[1],"-year = ", sprintf("%.3f", ROC.nomo$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",family='Times',
           x = 0.75, y = 0.25, size = 6,
           label = paste0("AUC at ",years[2],"-year = ", sprintf("%.3f", ROC.nomo$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",family='Times',
           x = 0.75, y = 0.15, size = 6,
           label = paste0("AUC at ",years[3],"-year = ", sprintf("%.3f", ROC.nomo$AUC[[3]])), color = "#E18727FF"
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
ggsave('14_ModeltimeROC.pdf',train_timeROC,w=7,h=7)
ggsave('14_ModeltimeROC.png',train_timeROC,w=7,h=7)


#DCA曲线
library(caret)
#BiocManager::install("ggDCA")
library(ggDCA)

variate <- c("age","weight","Stage","riskScore")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(variate,
                                         sep = '',
                                         collapse = '+')))
Nomogram <- rms::cph(cox_data_prog,
                     data = as.data.frame(df_all2))

d_train <- dca(Nomogram, times=100)

p1 <- ggplot(d_train)
p1
#dev.off()
ggsave("15_nomogram_dca.png", p1, width = 8, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("15_nomogram_dca.pdf", p1, width = 8, height = 6, units = "in", bg = "white")


model1<-coxph(Surv(OS.time,OS)~riskScore+age+weight+Stage,
              data=df_all2)

df1 <- ggDCA::dca(model1, 
                  times = c(12*1)
)

df2 <- ggDCA::dca(model1, 
                  times = c(12*2) # 1.5年，默认值是中位数
)

df3 <- ggDCA::dca(model1, 
                  times = c(12*3) # 1.5年，默认值是中位数
)

p1 <- ggplot(df1,       
             model.names="1-year",
             linetype =F, #线型
             lwd = 1.2)+   #线粗
  labs(title = '1-year')
p2 <- ggplot(df2,       
             model.names="2-year",
             linetype =F, #线型
             lwd = 1.2)+   #线粗
  labs(title = '2-year')

p3 <- ggplot(df3,       
             model.names="3-year",
             linetype =F, #线型
             lwd = 1.2)+   #线粗
  labs(title = '3-year')

library(patchwork)
p.cavd <- p1+p2+p3+
  plot_layout(ncol = 3)
p.cavd


pdf("15_DCA.pdf", width = 14, height = 6)
print(p.cavd)
dev.off()
png("15_DCA.png", width = 14, height = 6,res = 600,units = 'in')
print(p.cavd)
dev.off()


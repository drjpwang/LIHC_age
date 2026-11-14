rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./15_TIDE")){
  dir.create("./15_TIDE")
}
setwd("./15_TIDE")

getwd()

#读取风险得分数据
risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)
risk$sample <- rownames(risk)
risk1 <- risk[,c(ncol(risk)-2,ncol(risk)-1,ncol(risk))]

#TIDE分析 使用count数据（有对照样本）
#读取所有的样本
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-counts.txt",sep="\t",header=T,row.names = 1,
                stringsAsFactors = F,quote="",)
range(data)
colnames(data)=substring(colnames(data),1,12)

# expr_c <- colnames(train_t)
# for (i in 1:length(expr_c)) {
#   expr_c[i] <- gsub("\\.", "-", expr_c[i])
# }
# colnames(train_t) <- expr_c
#标准化
train_t <- data
normalize_data <-  t(apply(train_t, 1, function(x)x-(mean(x))))
write.table(data.frame("ID"=rownames(normalize_data),normalize_data),file="TIDE_input.txt",sep = "\t",quote = F,row.names = F)
#zip::zip("./TIDE_input.zip", normalize_data)

tide <- read.table("LIHC_TIDE_result.csv", header = T, row.names = 1, sep = ",")
#提取肿瘤样本
tide_c <- rownames(tide)
for (i in 1:length(tide_c)) {
  #tide_c[i] <- substr(tide_c[i], 1, 15)
  tide_c[i] <- gsub("\\.", "-", tide_c[i])
}
rownames(tide) <- tide_c

tide <- tide[risk1$sample,]
tide$sample <- rownames(tide)
#比较高低风险组的TIDE差异
tide_d <- merge(tide, risk1, by = "sample")
rownames(tide_d) <- tide_d$sample
tide_d$risk <- as.factor(tide_d$risk)

library(ggpubr)
tide_plot <- ggplot(tide_d, aes(x = risk, y = TIDE, fill = risk)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.size = 2, 
               linewidth = 0.8, color = "black") +
  geom_jitter(shape = 21, position = position_jitter(0.15), 
              size = 2, color = "black", alpha = 0.6) +
  
  scale_fill_manual(values = c("#dc143c", "#6495ed"), name = "Risk") +
  
  # 关键修复：移除position参数，使用默认分组逻辑
  stat_compare_means(method = "wilcox.test",  # 无需重复指定group=risk（已由x=risk决定）
                     label = "p.signif", cex = 5,
                     vjust = -0.5, hjust = 0.5,
                     label.x = 1.5, #可用于调整星号的位置
                     size = 12) +  # 仅保留位置调整参数
  
  ylim(min(tide_d$TIDE) * 0.95, max(tide_d$TIDE) * 1.15) +
  
  labs(
    # title = "TIDE Score",
    x = "Risk Group", 
    y = "TIDE Score",
    fill = "Risk Group"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15)),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    legend.position = "top",
    panel.grid = element_blank()
  )
tide_plot
ggsave('02_train_Immune_TIDE_boxplot.pdf',tide_plot,w=8,h=6)
ggsave('02_train_Immune_TIDE_boxplot.png',tide_plot,w=8,h=6)

#TIDE评分与风险评分的相关性
p1<-ggplot(tide_d,aes(x=riskScore,y=TIDE))+
  geom_point(aes(x=riskScore,y=TIDE,color =risk), size=2,alpha=0.9)+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ 
  stat_cor(data=tide_d, method = "spearman")+
  scale_fill_manual(values =c("black")) +
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="Risk Score",y="TIDE Score")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1

ggsave('03_train_RiskScore_TIDE_correlation.pdf',p1,w=8,h=6)
ggsave('03_train_RiskScore_TIDE_correlation.png',p1,w=8,h=6)

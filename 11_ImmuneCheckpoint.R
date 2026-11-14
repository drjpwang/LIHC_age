rm(list = ls())
setwd("/data/LIHC")
if (! dir.exists("./11_ImmuneCheckpoint")){
  dir.create("./11_ImmuneCheckpoint")
}
setwd("./11_ImmuneCheckpoint")

getwd()

#读取风险得分数据
risk <- read.table("../03_Survival/11_all_train_riskScore.txt", header = T)
risk$sample <- rownames(risk)
risk1 <- risk[,c(ncol(risk)-1,ncol(risk))]

#读取所有的样本
data=read.table("../00.data/TCGA/TCGA-LIHC-mRNA-tpm_log_case.txt",sep="\t",header=T,stringsAsFactors = F,quote="",check.names = F)
colnames(data)=substring(colnames(data),1,12)

cnames <- risk1$sample
train_t <- data[,cnames]

#免疫检查点
checkpoint <- c("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1",
                "HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRBS","HLA-E",
                "HLA-F","HLA-G","MICB","LAG3","PDCD1","PDCD1LG2","CD274","CD8A","CTLA4","BTLA",
                "IDO1","ADORA2A","CD27","CD40","CD28","ICOS","ICOSLG","HAVCR2","BTN3A1","BTN3A2",
                "SLAMF7","ENTPD1","PRF1")

train_dat <- train_t[rownames(train_t) %in% checkpoint,]

train_dat <- t(train_dat) %>% as.data.frame()
train_dat$sample <- rownames(train_dat)

train_data <- merge(train_dat, risk1, by = "sample")
rownames(train_data) <- train_data$sample

#根据高低风险组对疾病样本进行分组
colnames(train_data)[ncol(train_data)] <- "group"

train.dat_g <- gather(data = train_data,key = gene,value = Expression,-c("sample","group"))
train_plot <- ggplot(train.dat_g, aes(x=gene, y=Expression, fill = group)) +
  scale_fill_manual(values = c("#dc143c","#6495ed")) + 
  #ylim(4.2,11.3)+
  geom_boxplot()+
  stat_compare_means(aes(group = group),method = 'wilcox.test',label = "p.signif",cex=4)+
  # annotate(geom = "text", x = violin_dat1$gene, y = 15, size = 5,
  # label =as.character( violin_dat1$P.Value.signif)) +
  xlab("Gene Name") + 
  ylab("Expression") + 
  labs(title = "Survival_Expression_Immune")+ 
  theme(panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA))+
  theme(axis.text.x=element_text(face = "bold", color="black",angle = 70,vjust = 1, hjust = 1 ))+#angle = 70
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size=14),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size=12))+
  theme(plot.title = element_text(size=16))
train_plot
ggsave('01_train_Immune_boxplot.pdf',train_plot,w=14,h=8)
ggsave('01_train_Immune_boxplot.png',train_plot,w=14,h=8)

rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./20_CellChat")){
  dir.create("./20_CellChat")
}
setwd("./20_CellChat")

getwd()


#细胞通讯
#BiocManager::install("cpplot")
library(CellChat)
library(Seurat)
#读取巨噬细胞的单细胞数据
load("../19_scRNA-seq_inferCNV/LIHC_sc_anno.rda")

scRNA <- sc_dat2

pacman::p_load(CellChat, Seurat,ggplot2,tidyverse)


### PART J-1:Normal   ---------------------------------------------------------
#读入seurat处理后的rds文件
data.input<-scRNA@assays$RNA$data # normalized data matrix
meta <- scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$type)

meta <- meta[meta$cell_type != "Malignant cells",]
cell.use <- rownames(meta)[meta$type == "Normal"]
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
table(meta$cell_type)
unique(meta$cell_type)

##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type",)

##增加其他meta信息
cellchat <- addMeta(cellchat, meta = meta)

# 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
cellchat <- setIdent(cellchat, ident.use = "cell_type")

levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

###  基于配受体分析的数据库
###  包括人和老鼠的，这里我们用人
CellChatDB <- CellChatDB.human
###  作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
showDatabaseCategory(CellChatDB)

###  看一下CellChatDB的基本结构
dplyr::glimpse(CellChatDB$interaction)
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use

### 表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
### 首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedGenes(cellchat)
### 然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
cellchat <- identifyOverExpressedInteractions(cellchat)

# ??CellChat::projectData
# library(CellChat)
# ls("package:CellChat")
# packageVersion("CellChat")

# project gene expression data onto PPI network (optional)
# cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat)
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)
#df.net<-na.omit(df.net)
write.csv(df.net,'01.df.net.Control.csv',row.names = F,quote = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

## 相互作用数目
# pdf('03.number_of_interactions.pdf',w=6,h=5,family='Times')
netVisual_heatmap(cellchat)
# dev.off()
# png('03.number_of_interactions.png',w=500,h=400,family='Times')
# netVisual_heatmap(cellchat)
# dev.off()
#> Do heatmap based on a merged object
# netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object

saveRDS(cellchat,file = 'Control.rds')
#### circle
pdf('01.net_number.Control.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

png('01.net_number.Control.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()
#
pdf('02.net_weight.Control.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

png('02.net_weight.Control.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)


#需要指定受体细胞和配体细胞
net_bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE)+
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        # plot.title = element_text(size=20),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
ggsave('03.buble.Control.pdf',net_bubble,w=18,h=12)

ggsave('03.buble.Control.png',net_bubble,w=18,h=12)

### PART J-2:AMI   ---------------------------------------------------------
#读入seurat处理后的rds文件
rm(list = ls())
load("../05_scRNA-seq/LIHC_sc_anno.rda")

scRNA <- sc_dat2

pacman::p_load(CellChat, Seurat,ggplot2,tidyverse)

data.input<-scRNA@assays$RNA$data # normalized data matrix
meta <-  scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$type)

cell.use <-  rownames(meta)[meta$type == "Tumor"]
data.input <-  data.input[, cell.use]
meta <-  meta[cell.use, ]
table(meta$cell_type)
unique(meta$cell_type)

##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type",)
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息9
cellchat <- setIdent(cellchat, ident.use = "cell_type") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用小鼠
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半

dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use
##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# project gene expression data onto PPI network (optional)
# cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat)
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)
#df.net<-na.omit(df.net)
write.csv(df.net,'04.df.net.LIHC.csv',row.names = F,quote = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

## 相互作用数目
# pdf('03.number_of_interactions.pdf',w=6,h=5,family='Times')
netVisual_heatmap(cellchat)
# dev.off()
# png('03.number_of_interactions.png',w=500,h=400,family='Times')
# netVisual_heatmap(cellchat)
# dev.off()
#> Do heatmap based on a merged object
# netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object

saveRDS(cellchat,file = 'LIHC.rds')
# cellchat <- readRDS(file = 'Tumor.rds')
#### circle
pdf('05.net_number.LIHC.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

png('05.net_number.LIHC.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
#
pdf('06.net_weight.LIHC.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('06.net_weight.LIHC.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
levels(cellchat@idents)


# show all the significant interactions (L-R pairs)
# cellchat <- readRDS('AMI.rds')
#需要指定受体细胞和配体细胞
cellchat@data.signaling

net_bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE)+
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        # plot.title = element_text(size=20),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))
ggsave('07.buble.LIHC.pdf',net_bubble,w=24,h=14)

ggsave('07.buble.LIHC.png',net_bubble,w=16,h=12)


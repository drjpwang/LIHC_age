rm(list = ls())
setwd("/data/LIHC")
if(! dir.exists("./06_Keycell")){
  dir.create("./06_Keycell")
}
setwd("./06_Keycell")

getwd()


load("../05_scRNA-seq/LIHC_sc_anno.rda")
table(Idents(sc_dat2))
head(sc_dat2@meta.data)
table(sc_dat2$type)

lasso_geneids <- readRDS("../03_Survival/08_lasso_geneids.rds")
feaplot <- FeaturePlot(sc_dat2,
                       features = lasso_geneids,repel = T, label = F,
                       slot = "data", cols = c("lightgrey","red"), ncol = 2, order = T)
ggsave("03_feature_featureplot.pdf", feaplot, width = 10, height = 10)
ggsave("03_feature_featureplot.png", feaplot, width = 10, height = 10)

vlplot <- VlnPlot(sc_dat2, features = lasso_geneids, group.by = 'cell_type', pt.size = 0)
ggsave("01_feature_vlnplot.pdf", vlplot, width = 10, height = 8)
ggsave("01_feature_vlnplot.png", vlplot, width = 10, height = 8)

dotplot <- DotPlot(sc_dat2, features = lasso_geneids)+
  coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','skyblue','pink','red')) #颜色
ggsave("02_feature_dotplot.pdf", dotplot, width = 10, height = 8)
ggsave("02_feature_dotplot.png", dotplot, width = 10, height = 8)


de_celltype <- c("T cells", "NKT cells", "Macrophage", "NK cells", "Endothelial",
                 "Hepatocytes", "Plasma cells", "Fibroblast", "B cells", "Epithelial", "DC")
my_comparisons <- list(c("Tumor","Normal"))



# hubgene <- c("TK2")
pacman::p_load(reshape2,tidyverse,Seurat,rstatix,ggpubr,gridExtra)

gene <- lasso_geneids
gene <-as.data.frame(gene)

som.scobj <- sc_dat2
head(som.scobj@meta.data)
# som.scobj <- JoinLayers(som.scobj)
# Idents(som.scobj) <- som.scobj$cell_type
# som.scobj@assays$RNA$data["MANDP2",]
j <- 1
all_plots <- list()

# pacman::p_load(reshape2,tidyverse,Seurat,rstatix,ggpubr,gridExtra)
i = 1
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-c(8,10)]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 14)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 10, h = 12, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 2
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-c(7,10)]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 14)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 14, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 3
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-10]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 18)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 18, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 4
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-10]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 18)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 18, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 5
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-c(7,10)]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 14)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 14, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 6
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-10]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 18)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 18, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 7
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-c(7,10)]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 14)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 14, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()



i = 8
vln.df <- data.frame(exp = som.scobj@assays$RNA$data[gene[i,],])
vln.df$sample <- rownames(vln.df)
vln.df <- reshape::melt(vln.df, id = "sample")
colnames(vln.df)[2:3] <- c("symbol", "exp")
vln.df <- filter_if(vln.df, is.numeric, all_vars((.) > 0))

anno <- som.scobj@meta.data[, c("type", "cell_type")]
anno$sample <- rownames(anno)
vln.df <- inner_join(vln.df, anno, by = "sample")
vln.df$symbol <- paste0(gene[i, ], ' Expression Level')

# 检查每个 celltype 和 group 的样本数量
sample_counts <- vln.df %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 过滤掉样本数量不足的组
vln.df_filtered <- vln.df %>%
  group_by(cell_type, type) %>%
  filter(n() >= 2) %>%
  ungroup()

# 再次检查样本数量
sample_counts <- vln.df_filtered %>%
  group_by(cell_type, type) %>%
  summarise(count = n(), .groups = 'drop')

print(sample_counts)  # 打印样本数量以供检查

# 进行 Wilcoxon 检验
stat_cell <- vln.df_filtered %>%
  group_by(cell_type) %>%
  filter(n_distinct(type) > 1) %>%  # 确保每个 celltype 至少有两个不同的组
  wilcox_test(exp ~ type) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

plotcelltype <- de_celltype[-10]
plots <- list()

for (n in 1:length(plotcelltype)) {
  type <- plotcelltype[n]
  plot.dat <- vln.df_filtered[vln.df_filtered$cell_type == type, ]
  
  if (nrow(plot.dat) > 0) {  # 确保有数据
    position <- c(max(plot.dat$exp) + 0.5, max(plot.dat$exp) + 0.8, max(plot.dat$exp) + 1.1)
    
    p <- ggplot(plot.dat, aes(x = type, y = exp, fill = type)) +
      geom_violin(trim = FALSE, color = 'black') +
      stat_boxplot(geom = 'errorbar', width = 0.1, position = position_dodge(0.9)) +
      geom_boxplot(width = 0.4, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
      scale_fill_manual(values = c("#ed556a", "#5cb3cc"), name = "Group") +
      labs(title = type, x = "Group", y = "Expression", size = 20) +  # 设置坐标标题
      geom_signif(comparisons = my_comparisons,
                  test = "wilcox.test",
                  map_signif_level = TRUE,
                  y_position = position + 0.9 ,
                  textsize = 3) +
      ggtitle(plot.dat$cell_type[1]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", face = "bold", size = 12),
            axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.position = "top",
            panel.grid = element_line(color = 'gray', size = 0.2),
            panel.border = element_rect(colour = 'black', fill = NA, size = 1)) +
      guides(fill = 'none')
    
    plots[[n]] <- p
    all_plots[[j]] <- p
    j <- j + 1
  }
}
title <- text_grob(gene[i, ], size = 16, face = "bold")

pdf(file = paste0("0", i, ".", gene[i, ], "_exp.pdf"), w = 12, h = 18)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

png(file = paste0("0", i, ".", gene[i, ], "_exp.png"), w = 12, h = 18, units = 'in', res = 600)
grid.arrange(grobs = plots, ncol = 3, top = title)
dev.off()

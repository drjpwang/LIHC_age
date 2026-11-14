rm(list = ls())
setwd("/data/LIHC/")
if(! dir.exists("./17_inferCNV")){
  dir.create("./17_inferCNV")
}
setwd("./17_inferCNV")

getwd()

library(Seurat)
#创建Infercnv对象
#先下载JAGS并直接安装
# BiocManager::install("rjags")
library(rjags)
# BiocManager::install("infercnv")
library(infercnv)
options("Seurat.object.assay.version" = "v3")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = "../16_Hepatocytes/expfile.txt",
                                    annotations_file = "../16_Hepatocytes/groupfile.txt",
                                    delim = "\t", 
                                    gene_order_file = "../16_Hepatocytes/genefile.txt",
                                    ref_group_names = c("Fibroblast(non-malignant)"))
#运行infercnv（和正常癌旁细胞进行对照）

infercnv_result = infercnv::run(infercnv_obj,
                                cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir="/data/nas1/zhoujie/LIHC/17_inferCNV/inferCNV",  # dir is auto-created for storing outputs
                                cluster_by_groups=F,   # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                                plot_steps=FALSE,
                                denoise=F,
                                HMM=F,  # 是否基于HMM预测CNV,True的话时间很久
                                output_format = "pdf"
)

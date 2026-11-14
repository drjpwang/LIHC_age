# LIHC Age-related Study

This repository contains the R scripts and processed single-cell RNA-sequencing data associated with our study of age-related molecular and cellular characteristics in liver hepatocellular carcinoma (LIHC).

## Repository contents

The repository includes scripts for the following analyses:

### Survival and external validation

- `03_Survival.R`
- `04_GSE14520_SurvivalVerify.R`
- `04_ICGC_SurvivalVerify.R`
- `08_Cox_Nomograph.R`

### Single-cell RNA-seq analysis

- `05_scRNA-seq.R`
- `05_scRNA-seqCellEnrichment.R`
- `05_scRNA-seq_CellEnrichment_Artery-1 ECs.R`
- `06_Keycell.R`
- `16_Hepatocytes.R`
- `21_Endothelial.R`
- `21_Endothelial_ssgsea.R`

### Functional enrichment analysis

- `14_GSEA_GSVA.R`
- `14_GSEA_GSVA_GSE14520.R`
- `14_GSEA_GSVA_ICGC.R`

### Immune-related analysis

- `09_ImmuneInfiltration.R`
- `11_ImmuneCheckpoint.R`
- `11_ImmuneInfunction.R`
- `15_TIDE.R`
- `28_ImmuneTherapy.R`

### Copy-number variation and cell communication

- `17_inferCNV.R`
- `18_inferCNV_result.R`
- `19_scRNA-seq_inferCNV.R`
- `20_CellChat.R`

### Regulatory networks, trajectory and deconvolution

- `23_Pseudotime_Endothelial.R`
- `27_Pyscenic.R`
- `28_BayesPrism.R`

## Data availability

The processed single-cell RNA-sequencing objects are deposited in
Zenodo because they exceed the GitHub file-size limit.

- `LIHC_sc_anno.rda`: annotated LIHC single-cell RNA-seq object
- `Endothelial_sc_anno.rda`: annotated endothelial-cell subset

Zenodo DOI:
[https://doi.org/10.5281/zenodo.24210484]

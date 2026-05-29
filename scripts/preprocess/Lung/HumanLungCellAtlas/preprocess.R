library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'lung'
project_name <- 'HumanLungCellAtlas'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




project_seuratobj <- readRDS("lung_HumanLungCellAtlas_integration.rds")
colnames(project_seuratobj@meta.data)
table(project_seuratobj@meta.data$ann_level_1)




epi_project_seuratobj <- subset(project_seuratobj, subset = ann_level_1 == 'Epithelial')
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = ann_level_1 %in% 'Immune')
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




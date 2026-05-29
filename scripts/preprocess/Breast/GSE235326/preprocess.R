library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'breast'
project_name <- 'HumanBreastCellAtlas'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




project_seuratobj <- readRDS("clean.rds")




epi_project_seuratobj <- subset(project_seuratobj, subset = broad_cell_type %in% c('basal','lumhr','lumsec'))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = broad_cell_type %in% c('bcells','myeloid','tcells'))
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




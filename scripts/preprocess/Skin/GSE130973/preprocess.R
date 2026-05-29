library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'skin'
project_name <- 'GSE130973'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




project_seuratobj <- readRDS("GSE130973_seurat_analysis_lyko.rds")
project_seuratobj = UpdateSeuratObject(object = project_seuratobj) # Update Seurat object created by Seurat v3
View(project_seuratobj@meta.data)
colnames(project_seuratobj@meta.data)
table(project_seuratobj@meta.data$subj)
table(project_seuratobj@meta.data$integrated_snn_res.0.4)
table(project_seuratobj@meta.data$age)


cluster_ids <- project_seuratobj@meta.data$integrated_snn_res.0.4


cell_type_annotation <- rep(NA, length(cluster_ids))


cell_type_annotation[cluster_ids %in% c(8, 10)] <- "Pericytes"
cell_type_annotation[cluster_ids == 4] <- "Vascular EC"
cell_type_annotation[cluster_ids == 12] <- "Lymphatic EC"
cell_type_annotation[cluster_ids %in% c(5, 7, 15)] <- "Keratinocytes"
cell_type_annotation[cluster_ids == 14] <- "Melanocytes"
cell_type_annotation[cluster_ids == 6] <- "T-cells"
cell_type_annotation[cluster_ids %in% c(1, 2, 3, 9)] <- "Fibroblasts"
cell_type_annotation[cluster_ids == 11] <- "Erythrocytes"
cell_type_annotation[cluster_ids %in% c(0, 16, 13)] <- "Macrophages/DC"


project_seuratobj@meta.data$cell_type <- cell_type_annotation

table(project_seuratobj@meta.data$cell_type)

project_seuratobj@meta.data$donor_id <- project_seuratobj@meta.data$subj




epi_project_seuratobj <- subset(project_seuratobj, subset = cell_type %in% c('Keratinocytes','Melanocytes'))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = cell_type %in% c('Macrophages/DC','T-cells'))
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




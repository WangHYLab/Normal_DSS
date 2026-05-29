library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'digestive_system'
project_name <- 'GSE115469'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




project_seuratobj <- readRDS("GSE181987_joint_fish_human.rds")
View(project_seuratobj@meta.data)
colnames(project_seuratobj@meta.data)
table(project_seuratobj@meta.data$orig.ident)
table(project_seuratobj@meta.data$group)
unique(project_seuratobj@meta.data[,c('orig.ident','group')])
table(project_seuratobj@meta.data$seurat_clusters)


project_seuratobj <- subset(project_seuratobj,subset = group=='Human')
project_seuratobj$donor_id <- project_seuratobj$orig.ident


marker_cluster <- read.csv('marker_cluster.csv',header = T,row.names = 1)

anno_df <- data.frame(seurat_clusters = 0:20,
                      cell_type = unique(marker_cluster$cluster))


project_seuratobj@meta.data$seurat_clusters <- as.character(project_seuratobj@meta.data$seurat_clusters)
anno_df$seurat_clusters <- as.character(anno_df$seurat_clusters)


project_seuratobj@meta.data$cell_type <- anno_df$cell_type[match(project_seuratobj@meta.data$seurat_clusters, anno_df$seurat_clusters)]
table(project_seuratobj@meta.data$cell_type)




epi_project_seuratobj <- subset(project_seuratobj, subset = cell_type %in% c(paste0('zhHep',1:7),'zhChol/Hep','zhChol1','zhChol2'))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = cell_type %in% c('zhInf mac','zhNK','zhNon-inf mac','zhPlasma','zhαβ T','zhγδ T'))
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




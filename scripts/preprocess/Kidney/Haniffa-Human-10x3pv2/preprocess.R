library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'kidney'
project_name <- 'Haniffa-Human-10x3pv2'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))


project_seuratobj <- readRDS('clean.rds')
head(colnames(project_seuratobj)) # Check cell IDs

metadata <- read.csv('metadata.csv',header = T)
View(metadata)
colnames(metadata)


rownames(metadata) <- metadata$DropletID


project_seuratobj <- AddMetaData(project_seuratobj,metadata = metadata)


table(project_seuratobj@meta.data$Label)
table(project_seuratobj@meta.data$Experiment) # Only adult samples: PapRCC, RCC1, RCC2, VHL_RCC, RCC3
table(project_seuratobj@meta.data$TissueDiseaseState) # All normal
table(project_seuratobj@meta.data$Organ)
table(project_seuratobj@meta.data$PatientDiseaseState)
table(project_seuratobj@meta.data$Compartment) # 19565 Indistinct cells, Normal_Epithelium_and_Vascular_without_PT not separated
table(project_seuratobj@meta.data$ClusterID)
table(project_seuratobj@meta.data$Organ) # Only kidney selected, original data may include Ureter
cluster_info <- read.csv('cluster_info.csv',header = T) # Contains cell type for each cluster
View(cluster_info)
colnames(cluster_info)
table(cluster_info$Cell_type1)
head(cluster_info[,c('Cluster_ID','Cell_type1')])

matched_cell_type <- cluster_info$Cell_type1[match(project_seuratobj@meta.data$ClusterID, cluster_info$Cluster_ID)]
project_seuratobj@meta.data$Cell_type1 <- matched_cell_type
unique(project_seuratobj@meta.data[, c('ClusterID', 'Cell_type1','Compartment')]) # View combined cell type info
table(project_seuratobj@meta.data$Cell_type1)




project_seuratobj <- subset(project_seuratobj, subset = Experiment %in% c('PapRCC','RCC1','RCC2','RCC3','VHL_RCC'))
project_seuratobj$donor_id <- project_seuratobj$Experiment


epi_project_seuratobj <- subset(project_seuratobj, subset = (Cell_type1 == 'Nephron_epithelium' | Compartment == "Normal_Proximal_Tubules"))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = (Compartment == 'Normal_Immune' & Cell_type1 != "Junk" & Cell_type1 != "Private"))
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




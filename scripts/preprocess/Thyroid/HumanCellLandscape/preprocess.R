library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'thyroid'
project_name <- 'HumanCellLandscape'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))
source('~/BRCA/sc_function.R')


metadata_all <- read.csv('HCL_Fig1_cell_Info.csv',header = T)


donor_ids <- c('AdultThyroid1','AdultThyroid2')

for (donor in donor_ids){
  
  # Read rds and inspect data
  donorObj <- readRDS(paste(tissue,project_name,donor,'raw.rds',sep = '_'))
  View(donorObj@meta.data)
  
  # Add metadata
  metadata <- metadata_all[metadata_all$batch==donor,] # Extract sample subset
  rownames(metadata) <- sapply(strsplit(metadata$cellnames,'\\.'),`[`,2) # Use second part of cellnames (split by .) as row names
  donorObj <- AddMetaData(donorObj,metadata)
  table(donorObj@meta.data$celltype)
 
  # Preprocessing (already processed)------
  
  # Subset and save as rds
  
  # Extract epithelial and immune subsets
  epi_donorObj <- subset(donorObj, subset = celltype == 'Thyroid follicular cell')
  saveRDS(epi_donorObj,paste(tissue,project_name,donor,'epi_raw.rds',sep = '_'))
  
  imm_donorObj <- subset(donorObj, subset = celltype %in% c('Antigen presenting cell (RPS high)','B cell','B cell (Plasmocyte)','Dendritic cell','M2 Macrophage','Macrophage','Mast cell','Monocyte','Neutrophil','Proliferating T cell','T cell'))
  saveRDS(imm_donorObj,paste(tissue,project_name,donor,'imm_raw.rds',sep = '_'))
  
}



message("250318，
/home/chenshuting/drug_toxicity/pan_tissue/Scripts/Tabula/run_scpharm_project_epi.R,
        Batch run Tabula epithelial/tissue data")

WORK_DIR <- '/home/chenshuting/drug_toxicity/pan_tissue/'

library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)
library(scPharm)

threshold.s = -1.435637
threshold.r = 1.531656

datasets_info <- read.csv(paste0(WORK_DIR,'Tables/sample_info_tissue.csv'))
source('/home/chenshuting/drug_toxicity/pan_tissue/run_scpharm_function.R')

meta_donor <- 'donor_id'

# Get all tissues
tissue_list <- unique(datasets_info$Tissue)

# Iterate through each tissue----
for (tissue in tissue_list) {
  message(paste("Processing:", tissue))
  
  # Iterate through all projects----
  for (project_name in datasets_info$Project[datasets_info$Tissue == tissue]) {
    message(paste("Processing:", tissue, project_name))
    DATA_DIR <- paste('/home/chenshuting/drug_toxicity/pan_tissue', tissue, project_name, 'raw_data/', sep = '/')
    setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue', tissue, project_name, sep = '/'))
    
    dse_files <- list.files(pattern = paste0("^",tissue,"_", project_name, "_(.*)_Dse\\.csv$")) # For immune data, replace tissue with Blood
    if (length(dse_files) == 0) # Continue only if no DSE files exist (no prior results)
    {
      # Construct file paths
      epi_path <- paste0(DATA_DIR, paste(tissue, project_name, "epi_raw.rds", sep = "_"))
      project_path <- paste0(DATA_DIR, paste(tissue, project_name, "raw.rds", sep = "_"))
      
      # Read Seurat object (prefer epi, then project)
      seurat_obj <- NULL
      if (file.exists(epi_path)) {
        seurat_obj <- readRDS(epi_path)
        # Run scPharm on the loaded dataset
        scPharm_on_project(seurat_obj, tissue, project_name, meta_donor)
      } else if (file.exists(project_path)) {
        seurat_obj <- readRDS(project_path)
        # Run scPharm on the loaded dataset
        scPharm_on_project(seurat_obj, tissue, project_name, meta_donor)
      } 
    }
    
    
  }
}


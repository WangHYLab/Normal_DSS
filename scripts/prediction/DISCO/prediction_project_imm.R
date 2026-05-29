message("250327，
/home/chenshuting/drug_toxicity/pan_tissue/Scripts/run_scpharm_project_imm.R,
        Batch run DISCO immune-tissue data,
        server interrupted on 250306, need to continue")

WORK_DIR <- '/home/chenshuting/drug_toxicity/pan_tissue/'

library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)
library(scPharm)

threshold.s = -1.589578
threshold.r = 1.648195

datasets_info <- read.csv('/home/chenshuting/drug_toxicity/data_resource/normal_scrna_data/disco/DISCO_datasets.csv')
source('/home/chenshuting/drug_toxicity/pan_tissue/run_scpharm_function.R')

# Get all tissues
tissue_list <- unique(datasets_info$Tissue)

# Iterate through each tissue----
for (tissue in tissue_list) {
  message(paste("Processing:", tissue))
  
  # Iterate through all projects----
  for (project_name in datasets_info$project_name[datasets_info$Tissue == tissue]) {
    message(paste("Processing:", tissue, project_name))
    DATA_DIR <- paste('/home/chenshuting/drug_toxicity/pan_tissue', tissue, project_name, 'raw_data/', sep = '/')
    setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue', tissue, project_name, sep = '/'))
    
    dse_files <- list.files(pattern = paste0("^Blood_", project_name, "_(.*)_Dse\\.csv$")) # For epithelial data, replace Blood with tissue
    if (length(dse_files) == 0) # Continue only if no DSE files exist (not run before)
      {
      # Construct file paths
      imm_path <- paste0(DATA_DIR, paste(tissue, project_name, "imm_raw.rds", sep = "_"))
      
      # Read Seurat object
      seurat_obj <- NULL
      if (file.exists(imm_path)) {
        seurat_obj <- readRDS(imm_path)
        
        # Run scPharm on the loaded dataset
        scPharm_on_project(seurat_obj, 'Blood', project_name, 'sample_id')
      } 
    }

   
  }
}


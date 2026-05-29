message("250523，
/home/chenshuting/drug_toxicity/pan_tissue/Scripts/GSE227136/run_scpharm_project_imm.R,
        Batch run GSE227136 immune data")

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

project_name <- 'GSE227136'
tissue <- 'Lung'

source('/home/chenshuting/drug_toxicity/pan_tissue/run_scpharm_function.R')
meta_donor <- 'donor_id'

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
    scPharm_on_project(seurat_obj, 'Blood', project_name, meta_donor)
  } 
} 

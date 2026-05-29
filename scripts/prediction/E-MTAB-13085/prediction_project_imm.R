message("250610，
/home/chenshuting/drug_toxicity/pan_tissue/Scripts/E-MTAB-13085/run_scpharm_project_imm.R,
        Batch run E-MTAB-13085 immune data")

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

datasets_info <- read.csv(paste0(WORK_DIR,'Tables/sample_info_imm.csv'))
project_name <- 'E-MTAB-13085'

source('/home/chenshuting/drug_toxicity/pan_tissue/run_scpharm_function.R')

meta_donor <- 'donor_id'


tissue_list <- unique(datasets_info$Tissue[datasets_info$Project==project_name])

# Iterate through each tissue----
for (tissue in tissue_list) {
  message(paste("Processing:", tissue))
  
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
}


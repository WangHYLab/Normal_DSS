

WORK_DIR <- '/home/chenshuting/drug_toxicity/null_distribution/'

library(scPharm)
library(tidyr)
library(mixtools)

gdscdata_tissue <- readRDS('/home/chenshuting/drug_toxicity/scpharm/Rdata/gdscdata_tissue.rds')
seurat_objects <- readRDS(paste0(WORK_DIR,'Rdata/seurat_objects_for_imm_null_distribution.rds'))

meta_data_list <- list()
for (key in names(seurat_objects)){
  # Extract tissue and project names
  tissue_project <- unlist(strsplit(key, "_", fixed = TRUE))
  tissue <- tissue_project[1]
  project <- tissue_project[2]
  seurat_obj <- seurat_objects[[key]]
  
  # Run scpharm
  seurat_obj <- scPharmIdentify_tissue(seurat_obj, 
                                       type = "cellline", 
                                       tissue = 'Blood', # Use blood signature for immune
                                       cancer = 'pan',
                                       bulkdata = bulkdata,
                                       gdscdata = gdscdata_tissue)
  meta.data <- seurat_obj@meta.data
  meta.data <- meta.data[, grep("scPharm_nes", colnames(meta.data), value = T)]
  meta.data.long <- gather(meta.data, key = "Drug", value = "NES")
  meta_data_list[[key]] <- meta.data.long
  
}
saveRDS(meta_data_list,
        paste0(WORK_DIR,
               'Rdata/NullDist_list_project_imm.rds'))
meta.data.long <- do.call(rbind, meta_data_list)

out.mix <- normalmixEM(meta.data.long$NES)
threshold.s <- out.mix$mu[1] - out.mix$sigma[1]
threshold.r <- out.mix$mu[2] + out.mix$sigma[2]
print(threshold.s)
print(threshold.r)

saveRDS(list(NullDist = meta.data.long$NES, threshold_r = threshold.r, threshold_s = threshold.s),
        paste0(WORK_DIR,
               'Rdata/NullDist_list_imm.rds'))



meta_data_list <- readRDS(paste0(WORK_DIR,
                                 'Rdata/NullDist_list_project_imm.rds'))
meta.data.long <- do.call(rbind, meta_data_list)

out.mix <- normalmixEM(meta.data.long$NES,k=2)
threshold.s <- out.mix$mu[1] - out.mix$sigma[1]
threshold.r <- out.mix$mu[2] + out.mix$sigma[2]
print(threshold.s)
print(threshold.r)

saveRDS(list(NullDist = meta.data.long$NES, threshold_r = threshold.r, threshold_s = threshold.s),
        paste0(WORK_DIR,
               'Rdata/NullDist_list_imm.rds'))

# Compute density for the null distribution
null_density <- density(NullDist_list$NullDist)

# Store density information as a data frame
null_density_df <- data.frame(x = null_density$x, y = null_density$y)

# Save the density data to a file
saveRDS(null_density_df, file = paste0(WORK_DIR,
                                       "Rdata/null_density_imm.rds"))






attach_meta_to_seurat <- function(seurat_obj, meta_data) {
  cell_names <- colnames(seurat_obj)
  num_cells <- length(cell_names)
  
  drugs <- unique(meta_data$Drug)
  
  for (drug in drugs) {
    nes_values <- meta_data$NES[meta_data$Drug == drug]
    
    if (length(nes_values) != num_cells) {
      stop(paste("Cell number mismatch for", drug))
    }
    
    seurat_obj@meta.data[[drug]] <- nes_values
  }
  
  return(seurat_obj)
}

for (obj_name in names(seurat_objects)) {
  seurat_objects[[obj_name]] <- attach_meta_to_seurat(
    seurat_obj = seurat_objects[[obj_name]],
    meta_data = meta_data_list[[obj_name]]
  )
}





threshold.s.immune = -1.589578
threshold.r.immune = 1.648195

label_drug_response <- function(seurat_obj, drug_names, sensitive_threshold, resistant_threshold) {
  for (drug in drug_names) {
    nes_values <- seurat_obj@meta.data[[drug]]
    
    label <- dplyr::case_when(
      nes_values < sensitive_threshold ~ "sensitive",
      nes_values > resistant_threshold ~ "resistant",
      TRUE ~ "other"
    )
    
    label_col <- gsub("nes", "label", drug)
    seurat_obj@meta.data[[label_col]] <- label
  }
  
  return(seurat_obj)
}


for (obj_name in names(seurat_objects)) {
  
  drug_names <- grep("^scPharm_nes_", colnames(seurat_objects[[obj_name]]@meta.data), value = TRUE)
  
  seurat_objects[[obj_name]] <- label_drug_response(
    seurat_obj = seurat_objects[[obj_name]],
    drug_names = drug_names,
    sensitive_threshold = threshold.s.immune,
    resistant_threshold = threshold.r.immune
  )
}

names(seurat_obj@meta.data)


dse_list <- list()

for (key in names(seurat_objects)){
  
  # Extract tissue and project names-----
  tissue_project <- unlist(strsplit(key, "_", fixed = TRUE))
  tissue <- tissue_project[1]
  project <- tissue_project[2]
  seurat_obj <- seurat_objects[[key]]

  
  # Calculate DSE
  seurat_obj@meta.data$cell.label <- 'adjacent'  # Assuming a static update
  Dse <- scPharmDse(seurat_obj)
  
  # Update Seurat objects with scPharm results
  seurat_objects[[key]] <- seurat_obj
  
  # Save DSE
  Dse_data <- data.frame(Drug=rownames(Dse),
                         Dse=Dse$Dse)
  dse_list[[key]] <- Dse_data
}

saveRDS(seurat_objects,paste0(WORK_DIR,'Rdata/seurat_objects_for_null_distribution_scpharmed_imm.rds'))
saveRDS(dse_list,
        paste0(WORK_DIR,
               'Rdata/NullDist_list_project_dse_imm.rds'))

dse.data.long <- do.call(rbind, dse_list)


saveRDS(dse.data.long,
        paste0(WORK_DIR,
               'Rdata/dse.data.long_imm.rds')) # Epithelial stores threshold, but not calculated here, p90 likely based on 90%

# Compute density for the null distribution
null_density <- density(dse.data.long$Dse)

# Store density information as a data frame
null_density_df <- data.frame(x = null_density$x, y = null_density$y)

# Save the density data to a file
saveRDS(null_density_df, file = paste0(WORK_DIR,
                                       "Rdata/null_density_dse_imm.rds"))

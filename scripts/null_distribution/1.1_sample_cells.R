# Sample cells from 3 datasets per tissue to build null distribution.
# Sample 3 datasets per tissue, 300 cells per dataset. Adjust if insufficient data, ensure 900 cells total per tissue

WORK_DIR <- '/home/chenshuting/drug_toxicity/null_distribution/'

datasets_info <- read.csv('/home/chenshuting/drug_toxicity/data_resource/normal_scrna_data/disco/DISCO_datasets.csv')

library(dplyr)
library(Seurat)

# Get all tissues
tissue_list <- unique(datasets_info$Tissue)

# Sample datasets
# Store qualified Seurat objects
seurat_objects <- list()

# Iterate through each tissue
for (tissue in tissue_list) {
  message(paste("Processing:", tissue))
  
  # Prescreen projects with sufficient cells
  valid_projects <- list()
  
  # Iterate through all projects
  for (project_name in datasets_info$project_name[datasets_info$Tissue == tissue]) {
    DATA_DIR <- paste('/home/chenshuting/drug_toxicity/pan_tissue', tissue, project_name, 'raw_data/', sep = '/')
    
    # Construct file paths
    epi_path <- paste0(DATA_DIR, paste(tissue, project_name, "epi_raw.rds", sep = "_"))
    project_path <- paste0(DATA_DIR, paste(tissue, project_name, "raw.rds", sep = "_"))
    
    # Read Seurat object (prefer epi, then project)
    seurat_obj <- NULL
    if (file.exists(epi_path)) {
      seurat_obj <- readRDS(epi_path)
    } else if (file.exists(project_path)) {
      seurat_obj <- readRDS(project_path)
    } 
    
    # Ensure Seurat object exists and has more than 300 cells
    if (!is.null(seurat_obj) && ncol(seurat_obj) >= 300) {
      valid_projects[[project_name]] <- seurat_obj
    }
  }
  
  # Skip if no valid projects found
  if (length(valid_projects) == 0) {
    message(paste("No valid projects found for:", tissue))
    next
  }
  
  # Randomly select 3 qualified projects
  selected_projects <- sample(valid_projects, min(3, length(valid_projects)))
  
  # Store in Seurat list
  for (project_name in names(selected_projects)) {
    key <- paste(tissue, project_name, sep = "_")
    seurat_objects[[key]] <- selected_projects[[project_name]]
  }
}


print(names(seurat_objects))






cell_counts <- data.frame(Tissue = character(), Project = character(), Num_Cells = integer(), stringsAsFactors = FALSE)


for (key in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[key]]
  num_cells <- ncol(seurat_obj)  # Get cell count
  
  # Split Tissue and Project names
  tissue_project <- strsplit(key, "_")[[1]]
  tissue <- tissue_project[1]
  project_name <- paste(tissue_project[-1], collapse = "_")  # Handle multiple underscore-separated project names
  
  # Add to dataframe
  cell_counts <- rbind(cell_counts, data.frame(Tissue = tissue, Project = project_name, Num_Cells = num_cells))
}


cell_counts <- cell_counts[order(cell_counts$Tissue, cell_counts$Project), ]


print(cell_counts)



tissue <- 'Bladder'
project_names <- c('GSE134355')
for (project_name in project_names){
  DATA_DIR <- paste('/home/chenshuting/drug_toxicity/pan_tissue', tissue, project_name, 'raw_data/', sep = '/')
  
  # Construct file paths
  epi_path <- paste0(DATA_DIR, paste(tissue, project_name, "epi_raw.rds", sep = "_"))
  project_path <- paste0(DATA_DIR, paste(tissue, project_name, "raw.rds", sep = "_"))
  
  # Read Seurat object (prefer epi, then project)
  seurat_obj <- NULL
  if (file.exists(epi_path)) {
    seurat_obj <- readRDS(epi_path)
  } else if (file.exists(project_path)) {
    seurat_obj <- readRDS(project_path)
  } 
  key <- paste(tissue, project_name, sep = "_")
  seurat_objects[[key]] <- seurat_obj
}




cell_counts <- cell_counts %>%
  mutate(Sample_num_cells = ifelse(as.numeric(Num_Cells) < 300, as.numeric(Num_Cells), 300))
cell_counts$Sample_num_cells[cell_counts$Project == 'PRJNA662018' & cell_counts$Tissue == 'Bladder'] <- 808
cell_counts$Sample_num_cells[cell_counts$Project == 'GSE134520' & cell_counts$Tissue == 'Stomach'] <- 577
cell_counts$Sample_num_cells[cell_counts$Project == 'GSE159929' & cell_counts$Tissue == 'Stomach'] <- 323
write.csv(cell_counts,
          paste0(WORK_DIR,'Tables/sample_disco_cells_for_null_distribution.csv'),
          row.names = F)


for (key in names(seurat_objects)) {
  # Extract tissue and project names
  tissue_project <- unlist(strsplit(key, "_", fixed = TRUE))
  tissue <- tissue_project[1]
  project <- tissue_project[2]
  
  # Get corresponding cell count
  target_cells <- cell_counts %>%
    filter(Tissue == tissue & Project == project) %>%
    pull(Sample_num_cells)
  
  # Skip if no match found
  if (length(target_cells) == 0) {
    message(paste("No sampling info for:", key))
    next
  }
  
  # Perform cell sampling
  seurat_obj <- seurat_objects[[key]]
  if (ncol(seurat_obj) > target_cells) {
    sampled_cells <- sample(colnames(seurat_obj), target_cells)
    seurat_objects[[key]] <- subset(seurat_obj, cells = sampled_cells)
    message(paste("Sampled", target_cells, "cells from", key))
  } else {
    message(paste("Skipping", key, "as it has", ncol(seurat_obj), "cells"))
  }
}


sapply(seurat_objects, ncol)
saveRDS(seurat_objects,
        paste0(WORK_DIR,'Rdata/seurat_objects_for_null_distribution.rds'))


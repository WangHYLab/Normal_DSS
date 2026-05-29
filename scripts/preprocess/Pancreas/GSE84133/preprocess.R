library(dplyr)
library(Seurat)
library(Matrix)

tissue = 'pancreas'
project_name = 'GSE84133'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))

# Specify donor_ids and iterate
donor_ids <- c('human1','human2','human3','human4')

for (donor in donor_ids){
  
  # Read expression matrix and build Seurat object
  expr_matrix_file <- grep(paste0(donor,'_umifm_counts.csv.gz$'),dir(),value = TRUE) # value = TRUE ensures returning matched filenames instead of indices
  expr_matrix <- read.table(gzfile(expr_matrix_file), header = T, sep = ",")
  expr_matrix[1:5,1:5]
  
  # Process data (based on data format)
  expr_matrix_t <- t(expr_matrix[,c(-1,-2,-3)])
  colnames(expr_matrix_t) <- expr_matrix$X
  expr_matrix_sparse <- as(Matrix(as.matrix(expr_matrix_t), sparse = TRUE), "dgCMatrix") # Convert to sparse matrix
  donorObj <- CreateSeuratObject(counts = expr_matrix_sparse, project = project_name, min.cells = 3, min.features = 200) # Only 108 cells
  
  # Read annotation information and add metadata
  anno <- expr_matrix[,c(1,3)]
  rownames(anno) <- anno$X
  donorObj <- AddMetaData(donorObj,metadata = anno)
  
  # Plot QC metrics distribution to determine thresholds, remove outliers only
  donorObj[["percent.mt"]] <- PercentageFeatureSet(donorObj, pattern = "^MT-")
  donorObj[["percent.ribo"]] <- PercentageFeatureSet(donorObj, pattern = "^RP[SL]")
  VlnPlot(donorObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 2)
  
  # QC filtering
  nFeature_RNA_min = 200
  nFeature_RNA_max = 5000
  nCount_RNA_min = 400
  nCount_RNA_max = 20000
  percent.mt_max = 20
  percent.ribo_max = 50
  donorObj <- subset(donorObj, subset = nFeature_RNA > nFeature_RNA_min & nCount_RNA > nCount_RNA_min &
                       nFeature_RNA < nFeature_RNA_max & nCount_RNA < nCount_RNA_max &
                       percent.mt < percent.mt_max & percent.ribo < percent.ribo_max)
  
  # Preprocessing
  donorObj <- NormalizeData(donorObj)
  donorObj <- FindVariableFeatures(donorObj, selection.method = "vst", nfeatures = 3000)
  donorObj <- ScaleData(donorObj, verbose = FALSE)
  donorObj <- RunPCA(donorObj, npcs = 50, verbose = FALSE)
  
  # Separate epithelial and immune
  table(donorObj@meta.data$assigned_cluster)
  epi_donorObj <- subset(donorObj, subset = assigned_cluster %in% c('acinar','ductal'))
  saveRDS(donorObj,paste(tissue, project_name, donor, "epi_raw.rds", sep = '_'))
  
  imm_donorObj <- subset(donorObj, subset = assigned_cluster %in% c('macrophage','mast','t_cell'))
  saveRDS(donorObj,paste(tissue, project_name, donor, "imm_raw.rds", sep = '_'))
  print('raw.rds saved')
}



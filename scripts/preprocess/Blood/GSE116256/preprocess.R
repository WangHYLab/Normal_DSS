library(dplyr)
library(Seurat)
library(Matrix)

tissue = 'blood'
project_name = 'GSE116256'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))

# Specify donor_ids and iterate
donor_ids <- c('BM1','BM2','BM3','BM4','BM5-34p','BM5-34p38n')

for (donor in donor_ids){
  
  # Read expression matrix and build Seurat object
  expr_matrix_file <- grep(paste0(donor,'.dem.txt.gz$'),dir(),value = TRUE) # value = TRUE ensures returning matched filenames instead of indices
  expr_matrix <- read.table(gzfile(expr_matrix_file), header = T, sep = "\t")
  #View(expr_matrix)
  rownames(expr_matrix) <- expr_matrix$Gene
  expr_matrix$Gene <- NULL
  expr_matrix_sparse <- as(Matrix(as.matrix(expr_matrix), sparse = TRUE), "dgCMatrix") # Convert to sparse matrix
  donorObj <- CreateSeuratObject(counts = expr_matrix_sparse, project = project_name, min.cells = 3, min.features = 200) # Only 108 cells
  
  # Read annotation information and add metadata
  anno_file <- grep(paste0(donor,'.anno.txt.gz$'),dir(),value = TRUE)
  anno <- read.table(gzfile(anno_file), header = T, sep = "\t")
  #View(anno)
  rownames(anno) <- anno$Cell
  donorObj <- AddMetaData(donorObj,metadata = anno)
  
  # Plot QC metrics distribution to determine thresholds, remove outliers only
  donorObj[["percent.mt"]] <- PercentageFeatureSet(donorObj, pattern = "^MT-")
  donorObj[["percent.ribo"]] <- PercentageFeatureSet(donorObj, pattern = "^RP[SL]")
  VlnPlot(donorObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 2)
  
  # Annotation already done, skip below
  # # QC filtering
  # nFeature_RNA_min = 200
  # nFeature_RNA_max = 5000
  # nCount_RNA_min = 400
  # nCount_RNA_max = 20000
  # percent.mt_max = 20
  # percent.ribo_max = 50
  # donorObj <- subset(donorObj, subset = nFeature_RNA > nFeature_RNA_min & nCount_RNA > nCount_RNA_min &
  #                      nFeature_RNA < nFeature_RNA_max & nCount_RNA < nCount_RNA_max &
  #                      percent.mt < percent.mt_max & percent.ribo < percent.ribo_max)
  
  # Preprocessing
  donorObj <- NormalizeData(donorObj)
  donorObj <- FindVariableFeatures(donorObj, selection.method = "vst", nfeatures = 3000)
  donorObj <- ScaleData(donorObj, verbose = FALSE)
  donorObj <- RunPCA(donorObj, npcs = 50, verbose = FALSE)
  saveRDS(donorObj,paste(tissue, project_name, donor, "raw.rds", sep = '_'))
  print('raw.rds saved')
}



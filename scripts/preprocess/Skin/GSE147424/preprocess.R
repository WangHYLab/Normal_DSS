library(dplyr)
library(Seurat)

tissue = 'skin'
project_name = 'GSE147424'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))


gsm_ids <- dir('gsapub/ftp/pub/gen/GEND000109/GENDX000109/')
for (donor in gsm_ids){
  print(donor)
  
  # Build Seurat object
  raw_seurat_data <- Read10X(data.dir=paste0('gsapub/ftp/pub/gen/GEND000109/GENDX000109/',donor)) # Default is gene.column=1
  donorObj <- CreateSeuratObject(counts = raw_seurat_data, project = project_name, min.cells = 3, min.features = 200)
  donorObj@meta.data$donor_id <- donor
  
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
  saveRDS(donorObj,paste(tissue, project_name, donor, "raw.rds", sep = '_'))
  print('raw.rds saved')
}


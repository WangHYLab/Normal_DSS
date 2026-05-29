library(dplyr)
library(Seurat)

tissue = 'urogenital_system'
project_name = 'scAnalysisOfEndometriosis'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))



donor_ids <- c('GSM6102532_C01','GSM6102533_C02','GSM6102534_C03')
for (donor in donor_ids){
  print(donor)
  # Read Seurat object
  donorObj <- readRDS(paste(tissue,project_name,donor,'raw.rds',sep = '_'))
  colnames(donorObj@meta.data)
  
  # # Plot QC metrics distribution to determine thresholds, remove outliers only
  # donorObj[["percent.mt"]] <- PercentageFeatureSet(donorObj, pattern = "^MT-")
  # donorObj[["percent.ribo"]] <- PercentageFeatureSet(donorObj, pattern = "^RP[SL]")
  # VlnPlot(donorObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 2)
  # 
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
  # 
  # # Preprocessing
  # donorObj <- NormalizeData(donorObj)
  # donorObj <- FindVariableFeatures(donorObj, selection.method = "vst", nfeatures = 3000)
  # donorObj <- ScaleData(donorObj, verbose = FALSE)
  # donorObj <- RunPCA(donorObj, npcs = 50, verbose = FALSE)
  # saveRDS(donorObj,paste(tissue, project_name, donor, "raw.rds", sep = '_'))
}


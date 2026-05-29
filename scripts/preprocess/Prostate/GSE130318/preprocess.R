library(dplyr)
library(Seurat)

tissue = 'urogenital_system'
project_name = 'GSE130318'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))



# 

# files <- dir()
# 

# gsm_files <- grep('^GSM', files, value = TRUE)
# 


# 

# for (gsm_id in gsm_ids) {

#   dir.create(gsm_id)
#   

#   gsm_related_files <- gsm_files[grep(gsm_id, gsm_files)]
#   

#   for (file in gsm_related_files) {
#     if (grepl("barcodes", file)) {
#       file.rename(file, file.path(gsm_id, "barcodes.tsv.gz"))
#     } else if (grepl("features", file)) {
#       file.rename(file, file.path(gsm_id, "features.tsv.gz"))
#     } else if (grepl("matrix", file)) {
#       file.rename(file, file.path(gsm_id, "matrix.mtx.gz"))
#     }
#   }
# }
# 




file_path <- "GENS00014826/features.tsv.gz"


features_data <- read.table(gzfile(file_path), header = FALSE, sep = "\t", nrows = 5)


gsm_ids <- c('GENS00014827','GENS00014826')
for (donor in gsm_ids){
  print(donor)
  
  # Remove GRCh38_ prefix
  file_path <- paste0(donor,"/features.tsv.gz")
  features_data <- read.table(gzfile(file_path), header = FALSE, sep = "\t")
  features_data$V2 <- gsub("GRCh38_", "", features_data$V2) # Remove 'GRCh38_' prefix from column 2
  write.table(features_data, gzfile(file_path), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) # Save again
  
  # Build Seurat object
  raw_seurat_data <- Read10X(data.dir=donor,gene.column=2) # Default is gene.column=1
  donorObj <- CreateSeuratObject(counts = raw_seurat_data, project = project_name, min.cells = 3, min.features = 200)
  donorObj@meta.data$donor_id <- donor
  
  # Plot QC metrics distribution to determine thresholds, remove outliers only
  donorObj[["percent.mt"]] <- PercentageFeatureSet(donorObj, pattern = "^MT-")
  donorObj[["percent.ribo"]] <- PercentageFeatureSet(donorObj, pattern = "^RP[SL]")
  VlnPlot(donorObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 2)
  
  # QC filtering
  nFeature_RNA_min = 200
  nFeature_RNA_max = 6000
  nCount_RNA_min = 400
  nCount_RNA_max = 40000
  percent.mt_max = 50
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


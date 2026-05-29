library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'blood'
project_name <- 'SCP548'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




expr_matrix <- read.table(gzfile("scp_gex_matrix.csv.gz"), sep = ",", header = TRUE, row.names = 1)


expr_matrix[1:5,1:5]


expr_matrix_sparse <- as(Matrix(as.matrix(expr_matrix), sparse = TRUE), "dgCMatrix")


project_seuratobj <- CreateSeuratObject(counts = expr_matrix_sparse, project = project_name, min.cells = 3, min.features = 200)


metadata <- read.csv('scp_meta_updated.txt', sep = "\t", header = TRUE, row.names = 1)
#View(metadata)
colnames(metadata) # Metadata shows mitochondrial ratio and cell type, likely preprocessed
metadata <- metadata[-1,] # Remove first row (type)
table(metadata$Cell_Type)
table(metadata$donor_id)
table(metadata$disease__ontology_label) # All normal
head(rownames(metadata))
rownames(metadata) <- gsub("-", ".", rownames(metadata)) # Replace - with . to match project_seuratobj


project_seuratobj <- AddMetaData(project_seuratobj, metadata = metadata)

# Preprocessing
# Plot QC metrics distribution to determine thresholds, remove outliers only
project_seuratobj[["percent.mt"]] <- PercentageFeatureSet(project_seuratobj, pattern = "^MT-")
project_seuratobj[["percent.ribo"]] <- PercentageFeatureSet(project_seuratobj, pattern = "^RP[SL]")
VlnPlot(project_seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 2)

# # QC filtering
# nFeature_RNA_min = 200
# nFeature_RNA_max = 5000
# nCount_RNA_min = 400
# nCount_RNA_max = 20000
# percent.mt_max = 20
# percent.ribo_max = 50
# project_seuratobj <- subset(project_seuratobj, subset = nFeature_RNA > nFeature_RNA_min & nCount_RNA > nCount_RNA_min &
#                      nFeature_RNA < nFeature_RNA_max & nCount_RNA < nCount_RNA_max &
#                      percent.mt < percent.mt_max & percent.ribo < percent.ribo_max)

# Preprocessing
project_seuratobj <- NormalizeData(project_seuratobj)
project_seuratobj <- FindVariableFeatures(project_seuratobj, selection.method = "vst", nfeatures = 3000)
project_seuratobj <- ScaleData(project_seuratobj, verbose = FALSE)
project_seuratobj <- RunPCA(project_seuratobj, npcs = 50, verbose = FALSE)

saveRDS(project_seuratobj,paste(tissue,project_name,'raw.rds',sep = '_'))


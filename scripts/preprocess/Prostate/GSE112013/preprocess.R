library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'urogenital_system'
project_name <- 'GSE112013'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




expr_matrix <- read.table(gzfile("GSE112013_Combined_UMI_table.txt.gz"), sep = "\t", header = TRUE, row.names = 1)


expr_matrix[1:5,1:5]


expr_matrix_sparse <- as(Matrix(as.matrix(expr_matrix), sparse = TRUE), "dgCMatrix")


project_seuratobj <- CreateSeuratObject(counts = expr_matrix_sparse, project = project_name, min.cells = 3, min.features = 200)


metadata <- read.table('meta.txt', sep = " ", header = TRUE, row.names = 1)
View(metadata)
colnames(metadata)
rownames(metadata) <- gsub('-','.',rownames(metadata))
metadata$donor_id <- sapply(strsplit(rownames(metadata), "\\."), `[`, 1)


project_seuratobj <- AddMetaData(project_seuratobj, metadata = metadata)
#project_seuratobj$donor_id <- project_seuratobj$Sample



cluster_ids <- project_seuratobj@meta.data$Final_clusters


cell_type_annotation <- rep(NA, length(cluster_ids))


cell_type_annotation[cluster_ids == 1] <- "SSCs"
cell_type_annotation[cluster_ids == 2] <- "Differentiating S'gonia"
cell_type_annotation[cluster_ids == 3] <- "Early primary S'cytes"
cell_type_annotation[cluster_ids == 4] <- "Late primary S'cytes"
cell_type_annotation[cluster_ids == 5] <- "Round S'tids"
cell_type_annotation[cluster_ids == 6] <- "Elongated S'tids"
cell_type_annotation[cluster_ids %in% c(7,8)] <- "Sperm"
cell_type_annotation[cluster_ids == 9] <- "Macrophages"
cell_type_annotation[cluster_ids == 10] <- "Endothelial Cells"
cell_type_annotation[cluster_ids == 11] <- "Myoid Cells"
cell_type_annotation[cluster_ids == 12] <- "Sertoli Cells"
cell_type_annotation[cluster_ids == 13] <- "Leydig Cells"


project_seuratobj@meta.data$cell_type <- cell_type_annotation


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




epi_project_seuratobj <- subset(project_seuratobj, subset = Cell.Type %in% c('basal','Deuterosomal cells','Ionocytes','Multiciliated','Secretory','suprabasal'))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = Cell.Type %in% c('B cells','cDC','Mast cells','Monocytes','pDC','T cells'))
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




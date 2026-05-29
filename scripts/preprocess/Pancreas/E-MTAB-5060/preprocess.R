library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'pancreas'
project_name <- 'E-MTAB-5060'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




expr_matrix <- readMM("E-MTAB-5061.expression_tpm.mtx")
expr_matrix[1:5,1:5] # This is already a normalized and filtered matrix


cell_ids <- readLines("E-MTAB-5061.expression_tpm.mtx_cols")
head(cell_ids)


gene_ids <- read.table("E-MTAB-5061.expression_tpm.mtx_rows", header = FALSE, stringsAsFactors = FALSE)
gene_ids <- gene_ids$V1  # Extract first column as gene names


rownames(expr_matrix) <- gene_ids
colnames(expr_matrix) <- cell_ids

expr_matrix_sparse <- as(Matrix(as.matrix(expr_matrix), sparse = TRUE), "dgCMatrix")



project_seuratobj <- CreateSeuratObject(counts = expr_matrix_sparse, project = project_name, min.cells = 3, min.features = 200)


metadata <- read.table('E-MTAB-5061.sdrf.txt', sep = "\t", header = TRUE)
View(metadata)
colnames(metadata) # Metadata shows mitochondrial ratio and cell type, likely preprocessed
rownames(metadata) <- metadata$Comment.ENA_RUN.


project_seuratobj <- AddMetaData(project_seuratobj, metadata = metadata)
project_seuratobj$donor_id <- project_seuratobj$Characteristics..individual.


table(metadata$Characteristics..disease.)
project_seuratobj <- subset(project_seuratobj,subset = Characteristics..disease.=='normal')

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
project_seuratobj@assays$RNA$data <- project_seuratobj@assays$RNA$counts # Already a normalized and filtered matrix
#project_seuratobj <- NormalizeData(project_seuratobj)
project_seuratobj <- FindVariableFeatures(project_seuratobj, selection.method = "vst", nfeatures = 3000)
project_seuratobj <- ScaleData(project_seuratobj, verbose = FALSE)
project_seuratobj <- RunPCA(project_seuratobj, npcs = 50, verbose = FALSE)


table(project_seuratobj@meta.data$Comment..submitted.inferred.cell.type.)


epi_project_seuratobj <- subset(project_seuratobj, subset = Comment..submitted.inferred.cell.type. %in% c('acinar cell','ductal cell'))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))





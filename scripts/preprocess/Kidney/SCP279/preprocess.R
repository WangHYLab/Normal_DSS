library(scPharm)
library(dplyr)
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)
library(stringr)

tissue <- 'kidney'
project_name <- 'SCP279'
setwd(paste('/home/chenshuting/drug_toxicity/pan_tissue',tissue,project_name,'raw_data',sep = '/'))




expr_matrix <- read.table(gzfile("exprMatrix.tsv.gz"), sep = "\t", header = TRUE, row.names = 1)
expr_matrix_metro <- read.table(gzfile("exprMatrixSleMetro.tsv.gz"), sep = "\t", header = TRUE, row.names = 1)
expr_matrix_metro <- read.table(gzfile("exprMatrixSleMetro.tsv.gz"), sep = "\t", header = TRUE, row.names = 1)


expr_matrix[1:5,1:5]
expr_matrix_metro[1:5,1:5]


expr_matrix_sparse <- as(Matrix(as.matrix(expr_matrix), sparse = TRUE), "dgCMatrix")


project_seuratobj <- CreateSeuratObject(counts = expr_matrix_sparse, project = project_name, min.cells = 3, min.features = 200)


metadata <- read.csv('meta.tsv', sep = "\t", header = TRUE, row.names = 1)
View(metadata)
colnames(metadata) # Metadata shows mitochondrial ratio and cell type, likely preprocessed


project_seuratobj <- AddMetaData(project_seuratobj, metadata = metadata)
project_seuratobj$donor_id <- project_seuratobj$Sample




epi_project_seuratobj <- subset(project_seuratobj, subset = Cell.Type %in% c('basal','Deuterosomal cells','Ionocytes','Multiciliated','Secretory','suprabasal'))
saveRDS(epi_project_seuratobj,paste(tissue,project_name,'epi_raw.rds',sep = '_'))

imm_project_seuratobj <- subset(project_seuratobj, subset = Cell.Type %in% c('B cells','cDC','Mast cells','Monocytes','pDC','T cells'))
saveRDS(imm_project_seuratobj,paste(tissue,project_name,'imm_raw.rds',sep = '_'))




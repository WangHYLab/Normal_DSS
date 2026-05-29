


WORK_DIR <- '/home/chenshuting/drug_toxicity/null_distribution/'

library(scPharm)
library(tidyr)
library(mixtools)
library(ggplot2)

gdscdata_tissue <- readRDS('/home/chenshuting/drug_toxicity/scpharm/Rdata/gdscdata_tissue.rds')
seurat_objects <- readRDS(paste0(WORK_DIR,'Rdata/seurat_objects_for_null_distribution.rds'))

meta_data_list <- list()
dse_list <- list()

for (key in names(seurat_objects)){
  
  # Extract tissue and project names-----
  tissue_project <- unlist(strsplit(key, "_", fixed = TRUE))
  tissue <- tissue_project[1]
  project <- tissue_project[2]
  seurat_obj <- seurat_objects[[key]]
  
  # Run scPharm to calculate NES and DSE
  seurat_obj <- scPharmIdentify_tissue(seurat_obj, 
                                       type = "cellline", 
                                       tissue = tissue, 
                                       cancer = 'pan',
                                       bulkdata = bulkdata,
                                       gdscdata = gdscdata_tissue)
  seurat_obj@meta.data$cell.label <- 'adjacent'  # Assuming a static update
  Dse <- scPharmDse(seurat_obj)

  # Update Seurat objects with scPharm results
  seurat_objects[[key]] <- seurat_obj
  
  # Save NES
  meta.data <- seurat_obj@meta.data
  meta.data <- meta.data[, grep("scPharm_nes", colnames(meta.data), value = T)]
  meta.data.long <- gather(meta.data, key = "Drug", value = "NES")
  meta_data_list[[key]] <- meta.data.long
  
  # Save DSE
  Dse_data <- data.frame(Drug=rownames(Dse),
                         Dse=Dse$Dse)
  dse_list[[key]] <- Dse_data
}

saveRDS(seurat_objects,paste0(WORK_DIR,'Rdata/seurat_objects_for_null_distribution_scpharmed.rds'))
saveRDS(meta_data_list,
        paste0(WORK_DIR,
               'Rdata/NullDist_list_project.rds'))
saveRDS(dse_list,
        paste0(WORK_DIR,
               'Rdata/NullDist_list_project_dse.rds'))


meta.data.long <- do.call(rbind, meta_data_list)

out.mix <- normalmixEM(meta.data.long$NES,k=2)
threshold.s <- out.mix$mu[1] - out.mix$sigma[1]
threshold.r <- out.mix$mu[2] + out.mix$sigma[2]
print(threshold.s)
print(threshold.r)

saveRDS(list(NullDist = meta.data.long$NES, threshold_r = threshold.r, threshold_s = threshold.s),
        paste0(WORK_DIR,
               'Rdata/NullDist_list.rds'))

# Compute density for the null distribution
null_density <- density(NullDist_list$NullDist)

# Store density information as a data frame
null_density_df <- data.frame(x = null_density$x, y = null_density$y)

# Save the density data to a file
saveRDS(null_density_df, file = paste0(WORK_DIR,
                                       "Rdata/null_density.rds"))



dse_list <- readRDS(paste0(WORK_DIR,'Rdata/NullDist_list_project_dse.rds'))
dse.data.long <- do.call(rbind, dse_list)

sigma_hat <- sqrt(mean(dse.data.long$Dse^2)) # Assume Dse follows half-normal distribution

NullDist_list <- list(NullDist = dse.data.long$Dse, threshold = sigma_hat)
# draw_dsenulldist_density_plot(NullDist_list)
saveRDS(NullDist_list,
        paste0(WORK_DIR,
               'Rdata/NullDist_list_Dse.rds'))

# Compute density for the null distribution
null_density <- density(NullDist_list$NullDist)

# Store density information as a data frame
null_density_df <- data.frame(x = null_density$x, y = null_density$y)

# Save the density data to a file
saveRDS(null_density_df, file = paste0(WORK_DIR,
                                       "Rdata/null_density_dse.rds"))




# ggplot(dse.data.long, aes(x = Dse)) +
#   geom_density(fill = "lightblue", alpha = 0.5) +
#   labs(title = "Density Plot of Dse", x = "Dse", y = "Density") +
#   theme_minimal()
# 

# ggplot(dse.data.long, aes(x = log(Dse + 1))) +
#   geom_density(fill = "lightblue", alpha = 0.5) +
#   labs(title = "Density Plot of log(Dse + 1)", x = "log(Dse + 1)", y = "Density") +
#   theme_minimal()
# 

# p90 <- quantile(dse.data.long$Dse, 0.90)
# p95 <- quantile(dse.data.long$Dse, 0.95)
# p99 <- quantile(dse.data.long$Dse, 0.99)
# 

# ggplot(dse.data.long, aes(x = Dse)) +
#   geom_density(fill = "lightblue", alpha = 0.5) +

#   geom_vline(xintercept = p95, color = "red", linetype = "dashed") +

#   geom_vline(xintercept = p90, color = "blue", linetype = "dashed") +
#   labs(
#     title = "Density Plot of Dse",
#     x = "Dse",
#     y = "Density"
#   ) +
#   theme_minimal()

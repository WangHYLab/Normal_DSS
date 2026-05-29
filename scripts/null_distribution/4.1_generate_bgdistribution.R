

WORK_DIR <- '/home/chenshuting/drug_toxicity/null_distribution/'
DATA_DIR <- '/home/chenshuting/drug_toxicity/pan_tissue_landscape/'

library(scPharm)
library(tidyr)
library(mixtools)
library(ggplot2)



sample_tox_df <- readRDS(paste0(DATA_DIR,'Rdata/sample_tox_df/sample_tox_df.rds'))

# Compute density for the null distribution
null_density <- density(sample_tox_df$DSE)

# Store density information as a data frame
null_density_df <- data.frame(x = null_density$x, y = null_density$y)

saveRDS(null_density_df,file = paste0(WORK_DIR,"Rdata/bg_distribution/bg_density_dse.rds"))

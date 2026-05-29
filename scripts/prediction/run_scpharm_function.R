
# 

library(dplyr)
library(tidyselect)
library(Seurat)
library(scPharm)
library(Matrix)
library(pROC)
library(reshape2)
library(scales)
library(ggplot2)

gdscdata_tissue <- readRDS('/home/chenshuting/drug_toxicity/scpharm/Rdata/gdscdata_tissue.rds')


threshold.s.tissue = -1.435637
threshold.r.tissue = 1.531656


threshold.dse.tissue = 0.04483429


threshold.s.immune = -1.589578
threshold.r.immune = 1.648195



scPharm_on_donor <- function(donorObj, tissue, project_name, donor) {
  # Step 1: Update donorObj with scPharmIdentify
  donorObj <- scPharmIdentify_tissue(donorObj, 
                                       type = "cellline", 
                                       tissue = tissue, 
                                       cancer = 'pan',
                                       bulkdata = bulkdata,
                                       gdscdata = gdscdata_tissue,
                                       threshold.s = threshold.s,
                                       threshold.r = threshold.r)
  donorObj@meta.data$cell.label <- 'adjacent'  # Assuming a static update
  
  # Step 2: Compute and Save Dse
  Dse <- scPharmDse(donorObj)
  colnames(Dse) <- c('DRUG_ID', 'DRUG_NAME', paste0(donor, '_', project_name))
  write.csv(Dse,
          file = paste(tissue, project_name, donor, "Dse.csv", sep = '_'),
          row.names = TRUE)
  
  # Step 3: Save the seuratObj with scPharm label and other associated data
  saveRDS(donorObj, 
          file = paste(tissue, project_name, donor, "scpharm.rds", sep = '_'))
  
  # Step 4: Save scPharm meta data
  scPharm_meta <- donorObj@meta.data[, grepl("^scPharm", colnames(donorObj@meta.data)), drop = FALSE]
  saveRDS(scPharm_meta, 
          file = paste(tissue, project_name, donor, "scpharmMeta.rds", sep = '_'))
  
  
}
# scPharm_on_donor(donorObj, tissue, project_name, donor)


scPharm_on_project <- function(project_seuratobj, tissue, project_name, meta_donor) {
  # Split Seurat object
  donor.list <- SplitObject(project_seuratobj, split.by = meta_donor)
  
  for (donor in names(donor.list)) {
    donorObj <- donor.list[[donor]]
    
    tryCatch({
      # Check for intermediate donor result files
      donor_result_file <- paste(tissue, project_name, donor, "scpharm.rds", sep = '_')
      
      if (file.exists(donor_result_file)) {
        print(paste("existing results for donor:", donor))
      } else {
        print(paste("Processing donor:", donor))
        
        # Call scPharm_on_donor function
        scPharm_on_donor(donorObj, tissue, project_name, donor)
      }
      
    }, error = function(e) {
      # Code to execute if error occurs
      message(paste("Error processing donor:", donor, "; Skipping to next. Error: ", e$message))
    })
  }
}

# scPharm_on_project(project_seuratobj, tissue, project_name, meta_donor)


merge_dse_files <- function(tissue, project_names, if_bloodsig = FALSE) {
  
  # Define output filename
  if (if_bloodsig) {
    output_file <- paste(tissue, 'imm',"Dse.rds", sep = '_')
  } else {
    output_file <- paste(tissue, "Dse.rds", sep = '_')
  }
  
  # Initialize empty list to store data frames
  data_list <- list()
  
  # Iterate through each project_name to get _Dse.csv files
  for (project_name in project_names) {
    
    # Construct path for each project_name folder
    project_path <- file.path('/home/chenshuting/drug_toxicity/pan_tissue', 
                              tissue,project_name)
    
    # Set different filename patterns based on if_bloodsig parameter
    if (if_bloodsig) {
      pattern <- paste0("blood_", project_name, "_.*_Dse.csv")
    } else {
      pattern <- paste0(tissue, "_", project_name, "_.*_Dse.csv")
    }
    
    # Get matching _Dse.csv files in this path
    file_list <- list.files(path = project_path, pattern = pattern, full.names = TRUE)
    
    # Iterate through each file, read data and store in list
    for (file_path in file_list) {
      # Read CSV file
      df <- read.csv(file_path, header = TRUE)
      
      # Set row names and remove first column
      rownames(df) <- df[, 1]
      df <- df[, -1]
      
      # Add to data list
      data_list[[file_path]] <- df
    }
  }
  
  # Merge all data frames
  Dse_merged <- Reduce(function(x, y) {
    merge(x, y, by = c("DRUG_ID", "DRUG_NAME"), all = TRUE)
  }, data_list)
  
  # Convert DRUG_ID to string format
  Dse_merged$DRUG_ID <- as.character(Dse_merged$DRUG_ID)
  
  # Set row names
  rownames(Dse_merged) <- paste('scPharm_label', Dse_merged$DRUG_ID, Dse_merged$DRUG_NAME, sep = '_')
  
  # Save merged data frame as .rds file in current tissue directory
  output_path <- file.path('/home/chenshuting/drug_toxicity/pan_tissue', 
                           tissue, output_file)
  saveRDS(Dse_merged, file = output_path)
  
  print(paste("Merged Dse data saved as:", output_file))
}

# merge_dse_files(tissue, project_names,if_bloodsig = FALSE)


merge_meta_files <- function(tissue, project_names, if_bloodsig = FALSE) {
  
  setwd(file.path('/home/chenshuting/drug_toxicity/pan_tissue', tissue))
  
  # Define output filename, merge meta files from all projects, output to tissue level
  if (if_bloodsig) {
    output_file <- paste(tissue, 'imm',"scpharmMeta.rds", sep = '_')
  } else {
    output_file <- paste(tissue, "scpharmMeta.rds", sep = '_')
  }
  
  # Initialize empty list to store meta data frames
  meta_list <- list()
  
  # Iterate through each project_name
  for (project_name in project_names) {
    
    # Construct path for each project_name folder
    project_path <- file.path(getwd(), project_name)
    
    # Set different filename patterns based on if_bloodsig parameter
    if (if_bloodsig) {
      pattern <- paste0("Blood_", project_name, "_.*_scpharmMeta.rds")
    } else {
      pattern <- paste0(tissue, "_", project_name, "_.*_scpharmMeta.rds")
    }
    
    # Get matching _scpharmMeta.rds files in this path
    file_list <- list.files(path = project_path, pattern = pattern, full.names = TRUE)
    
    # Iterate through each file, read meta info and store in list
    for (file_path in file_list) {
      # Read Seurat object metadata
      meta_data <- readRDS(file_path)
      
      # Add tissue information
      meta_data$tissue <- tissue
      
      # Add project_name information
      meta_data$project_name <- project_name
      
      # Extract donor info and add to meta data frame
      donor <- sub(paste0("^", tissue, "_", project_name, "_(.*)_scpharmMeta.rds$"), "\\1", basename(file_path))
      meta_data$donor <- donor
      
      # Add to meta list
      meta_list[[file_path]] <- meta_data
    }
  }
  
  # Merge all meta data frames
  merged_meta <- bind_rows(meta_list)
  
  # Save merged meta data frame as .rds file
  saveRDS(merged_meta, file = output_file)
  
  print(paste("Merged meta data saved as:", output_file))
}


# merge_meta_files(tissue, project_names, if_bloodsig = FALSE)



library(ggplot2)

draw_nulldist_density_plots <- function(NullDist_list, title_text = "Non-immune null distribution") {
  # Extract data and thresholds
  NullDist <- NullDist_list$NullDist
  threshold_r <- NullDist_list$threshold_r
  threshold_s <- NullDist_list$threshold_s
  
  df <- data.frame(NullDist = NullDist)
  
  # Calculate maximum density for text annotation height
  max_dens <- max(density(df$NullDist)$y)
  
  ggplot(df, aes(x = NullDist)) +
    # Plot filled density curve
    geom_density(fill = "lightblue", alpha = 0.5, color = "blue", size = 0.8) +
    
    # Plot red dashed threshold lines
    geom_vline(xintercept = threshold_r, linetype = "dashed", color = "red", size = 0.8) +
    geom_vline(xintercept = threshold_s, linetype = "dashed", color = "red", size = 0.8) +
    
    # Optimize text annotation: place on both sides of dashed line to avoid edge cutoff
    annotate("text", 
             x = threshold_s - 0.2, 
             y = max_dens * 0.7, 
             label = paste0("Threshold_s = ", round(threshold_s, 3)), 
             hjust = 1, size = 5) + 
    annotate("text", 
             x = threshold_r + 0.2, 
             y = max_dens * 0.7, 
             label = paste0("Threshold_r = ", round(threshold_r, 3)), 
             hjust = 0, size = 5) + 
    
    # Set labels
    labs(title = title_text,
         x = "NES value",
         y = "Density") +
    
    # Style modification: clean academic style, no grid, with axes
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),   # Remove major grid
      panel.grid.minor = element_blank(),   # Remove minor grid
      panel.border = element_blank(),       # Remove border
      axis.line = element_line(color = "black"), # Show only X and Y axis lines
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Center and bold title
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 12, color = "black")
    )
}




draw_dsenulldist_density_plot <- function(NullDist_list) {
  # Plot Dse density plot, input is a list containing Dse values and a threshold
  NullDist <- NullDist_list$NullDist
  threshold <- NullDist_list$threshold  # Only one threshold
  
  df <- data.frame(Dse = NullDist)
  
  ggplot(df, aes(x = Dse)) +
    geom_density(fill = "lightblue", alpha = 0.5, color = "blue") +  # Plot density curve
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", size = 1) +  # Mark threshold
    annotate("text", 
             x = threshold + 0.02,  # Adjust position to avoid overlap
             y = max(density(df$Dse)$y) / 2, 
             label = paste0("Threshold = ", round(threshold, 3)), 
             hjust = 0, size = 5) +  # Right align
    labs(title = "Density Plot of Dse",
         x = "Dse Values",
         y = "Density") +
    theme_minimal()
}


draw_density_plots <- function(xlab, bg_values, target_values) { # Plot distribution density: grey background + red target distribution
  # Plot null distribution density
  plot(density(bg_values), col = "#999999", main = "",
       lwd = 1.5, lty = 2, cex.axis = 1.2, font.main = 1, tcl = -0.25, ylab = "", xlab = xlab,
       cex.lab = 1.5, bty = "l", las = 1, 
       ylim = c(0, max(max(density(bg_values)$y),max(density(target_values)$y))) # Y limit is max of both densities
  )
  
  # Add target_values density line
  lines(density(target_values), col = "#FF0000", lwd = 1.5)
}

draw_density_plots_bgdse <- function(title, bg_density_df, target_values, threshold = NULL) { 
  # Calculate target distribution density
  foreground_density <- density(target_values)
  foreground_density_df <- data.frame(x = foreground_density$x, y = foreground_density$y)
  
  p <- ggplot() +
    # Background density
    geom_line(data = bg_density_df, aes(x = x, y = y), color = "grey", size = 1, alpha = 0.5) +
    
    # Foreground density
    geom_line(data = foreground_density_df, aes(x = x, y = y), color = "red", size = 1, alpha = 0.8) +
    
    labs(title = title,
         x = "Dse Values",
         y = "Density") +
    theme_minimal()
  
  # Add threshold line if not empty
  if (!is.null(threshold)) {
    p <- p + geom_vline(xintercept = threshold, linetype = "dashed", color = "grey", size = 1)
  }
  
  return(p)
}

draw_density_plots_bgnes <- function(title, bg_density_df, target_values, threshold_r,threshold_s) { # Plot distribution density: stored grey background + red target values
  # Compute density for foreground distribution
  foreground_density <- density(target_values)
  foreground_density_df <- data.frame(x = foreground_density$x, y = foreground_density$y)
  
  # Overlay the foreground density on the stored background density
  ggplot() +
    # Background density (precomputed)
    geom_line(data = bg_density_df, aes(x = x, y = y), color = "grey", size = 1, alpha = 0.5) +
    
    # Foreground density
    geom_line(data = foreground_density_df, aes(x = x, y = y), color = "red", size = 1, alpha = 0.8) +
    
    # Thresholds
    geom_vline(xintercept = threshold_r, linetype = "dashed", color = "grey", size = 1) +
    geom_vline(xintercept = threshold_s, linetype = "dashed", color = "grey", size = 1) +
    
    # Labels and theme
    labs(title = title,
         x = "NES Values",
         y = "Density") +
    theme_minimal()
}


# # Compute density for the null distribution
# null_density <- density(NullDist_list$NullDist)
# 
# # Store density information as a data frame
# null_density_df <- data.frame(x = null_density$x, y = null_density$y)
# 
# # Save the density data to a file (optional)
# saveRDS(null_density_df, file = paste0(DATA_DIR,
#                                        "Rdata/null_density.rds"))

# bg_density_df <- readRDS('/home/chenshuting/drug_toxicity/null_distribution/Rdata/null_density.rds')

# tissue <- 'Bladder'
# target_values <- NullDist_df$NES[NullDist_df$Tissue==tissue]
# 

# title <- paste0('NES distribution of ', tissue)
# 
# draw_density_plots_bgnes(title,bg_density_df,target_values)



plot_stacked_bar <- function(seurat.combined, meta_var1, meta_var2) { # Plot stacked bar chart with combined conditions
  # Calculate count/proportion of meta_var2 within each meta_var1
  counts_table <- table(seurat.combined@meta.data[[meta_var1]], seurat.combined@meta.data[[meta_var2]])
  proportions_table <- prop.table(counts_table, margin = 1)
  
  # Create stacked bar chart
  p <- ggplot(data = as.data.frame(proportions_table), aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity") +
    labs(x = meta_var1, y = "proportion", fill = meta_var2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase x-axis text size and adjust angle
          axis.text.y = element_text(size = 12),  # Increase y-axis text size
          axis.title.x = element_text(size = 14),  # Increase x-axis title size
          axis.title.y = element_text(size = 14),  # Increase y-axis title size
          legend.title = element_text(size = 12),  # Increase legend title size
          legend.text = element_text(size = 12))  # Increase legend text size
  
  return(p)
}
# plot_stacked_bar(seurat.combined, meta_var1, meta_var2)

plot_roc_curve <- function(data, response_var, predictor_var='RRA_score') { # Evaluate RRA_score performance using ADR database as label
  # Dynamically select response and predictor variables
  response <- data[[response_var]]
  predictor <- data[[predictor_var]]
  
  # Check if response variable has sufficient levels
  if (length(unique(response)) < 2) {
    print(paste("Error: The response variable", response_var, "must have at least two levels."))
    return(NULL)
  }
  
  # Reverse predictor values (lower RRA_score indicates stronger ADR correlation)
  predictor <- -predictor
  
  # Calculate ROC object
  roc_obj <- roc(response, predictor)
  print(roc_obj)
  
  # Extract AUC value
  auc_value <- auc(roc_obj)
  print(paste("AUC:", auc_value))
  
  # Plot ROC curve
  roc_curve <- ggroc(roc_obj)
  
  # Visualize using ggplot2
  ggplot(roc_curve$data, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "darkorange", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "navy") +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title = response_var,
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC =", round(auc_value, 2)), color = "black", size = 5)
}

# plot_roc_curve(Dse_adrsubset, 'ADR_related1')

update_seurat_metadata <- function(seurat_obj, new_metadata, columns) { # Update metadata from another Seurat object, add if not exists, update if exists
  # Get Seurat object metadata
  seurat_metadata <- seurat_obj@meta.data
  
  # Iterate through each column name
  for (column in columns) {
    if (column %in% colnames(seurat_metadata)) {
      # If column exists, only update cells present in new metadata
      cells_to_update <- intersect(rownames(seurat_metadata),rownames(new_metadata))
      seurat_metadata[cells_to_update, column] <- new_metadata[cells_to_update, column]
    } else {
      # If column doesn't exist, add column and initialize with NA
      seurat_metadata[[column]] <- NA
      # Update values for cells present in new metadata
      cells_to_update <- intersect(rownames(seurat_metadata),rownames(new_metadata))
      seurat_metadata[cells_to_update, column] <- new_metadata[cells_to_update, column]
    }
  }
  
  # Update Seurat object metadata
  seurat_obj@meta.data <- seurat_metadata
  return(seurat_obj)
}

# seurat_obj <- update_seurat_metadata(seurat.combined, metadata_np, annocolumns)


visualize_drug_info <- function(drug_info_subset) {
  
  # Convert dataframe to long format
  drug_info_long <- reshape2::melt(drug_info_subset, id.vars = c("DRUG_ID", "DRUG_NAME"))
  
  # Visualization function: map 1 to red, 0 to green, NA to white, gradient in between
  value_to_color <- function(val) {
    if (is.na(val)) {
      return("white")
    } else if (val == 1) {
      return("red")
    } else if (val == 0) {
      return("green")
    } else {
      return(scales::col_numeric(palette = c("green", "red"), domain = c(0, 1))(val))
    }
  }
  
  # Map colors for each value
  drug_info_long$color <- sapply(drug_info_long$value, value_to_color)
  
  # Create and return ggplot chart
  p <- ggplot(drug_info_long, aes(x = variable, y = interaction(DRUG_ID, DRUG_NAME))) +
    geom_point(aes(fill = color), shape = 21, size = 6, color = "black") +  # Adjust circle size
    scale_fill_identity() +
    labs(x = "", y = "Drug ID and Name", title = "Drug Safety Metrics Visualization") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),  
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),  # Remove grid lines
      plot.background = element_rect(fill = "white", color = NA),  # Set plot background to white
      panel.background = element_rect(fill = "white", color = NA)  # Set panel background to white
    ) +
    scale_y_discrete(limits = rev(unique(interaction(drug_info_subset$DRUG_ID, drug_info_subset$DRUG_NAME))))
  
  return(p)
}

#p <- visualize_drug_info(drug_info_subset)
# ggsave(filename = paste(tissue, project_name, "drug_safety_highreso.png", sep = '_'), plot = p, width = 12, height = 80, dpi = 300, limitsize = FALSE)



extract_drug_info <- function(drug_markers_list, columns) {
  # Use lapply to iterate through list and extract specified columns
  drug_info <- lapply(drug_markers_list, function(drug) {
    # Set default values for each column to prevent missing columns
    data <- lapply(columns, function(col) if (!is.null(drug[[col]])) drug[[col]] else NA)
    names(data) <- columns
    return(as.data.frame(data, stringsAsFactors = FALSE))
  })

  # Convert list to dataframe
  drug_info_df <- do.call(rbind, drug_info)
  return(drug_info_df)
}


#drug_info_df <- extract_drug_info(drug_markers_list, c('drug_id','has_overlap_genes','has_overlap_pathways','has_overlap_offgenes','has_overlap_offpathways'))


create_newcol <- function(df, col1, col2, new_col_name) {
  df <- df %>%
    mutate(!!new_col_name := ifelse(is.na(.[[col1]]) | is.na(.[[col2]]), 
                                    NA, 
                                    ifelse(.[[col1]] & .[[col2]], TRUE, FALSE)))
  return(df)
}

# target_info_df <- create_newcol(target_info_df, 'has_overlap_genes', 'has_overlap_pathways', 'match_MOA')


RenameGenesSeurat <- function(obj ,newnames ) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. 
  # It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

# obj=RenameGenesSeurat(obj = sce, newnames = ids$SYMBOL)


#########################
# Load necessary libraries
#
#
library(dplyr)
library(tidyr)
library(randomForest)
library(mgcv)
library(ggplot2)
library(xgboost)
library(boot)
library(iml)
library(caret)
library(glmnet)  # For Elastic Net Regression
library(gbm)      # For Gradient Boosting Machines
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(msigdbr)


# Set working directory (update path as needed)
setwd("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP")



# Load / Read necessary files
pheno <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/ROSMAPphenometh.csv")
samples <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/samplesROSMAP.csv")
load("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/regressedROSMAP.rda")
probes <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/450kinfo.csv")
#
read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/GAM_results.csv")
read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/GAM_results_with_cv.csv")
#
##########
  




##########################
#STEP 1: ORGANISE THE DATA
#
#
  # Convert betas.dasen to data frame and assign row names as a new column IlmnID
  betas_df <- as.data.frame(data.reg.ROSMAP)
  betas_df$IlmnID <- rownames(data.reg.ROSMAP)

  # Ensure IlmnID in both betas_df and probes are character for matching
  betas_df$IlmnID <- as.character(betas_df$IlmnID)
  probes$IlmnID <- as.character(probes$IlmnID)

  # Select relevant columns in each dataset
  probes2 <- probes %>% select(IlmnID, UCSC_RefGene_Name, CHR, MAPINFO, UCSC_RefGene_Group)
  pheno2 <- pheno %>% select(projid, Full_Sentrix, braaksc)
  samples2 <- samples %>% select(Full_Sentrix)

  # Merge Phenometh and Samples using Full_Sentrix as the identifier
  phenotype_data <- merge(pheno2, samples2, by = "Full_Sentrix", all.x = TRUE)

  # Reshape betas_df from wide to long format for merging
  betas_long <- betas_df %>%
   pivot_longer(
     cols = -IlmnID,  # Keep IlmnID as is, pivot all other columns
     names_to = "Full_Sentrix",  # Convert column names (sample IDs) to "Full_Sentrix"
      values_to = "MethylationValue"  # Methylation values go to "MethylationValue"
    )

  # Merge with phenotype data on Full_Sentrix
  merged_data <- betas_long %>%
  inner_join(phenotype_data, by = "Full_Sentrix")

  # Convert to wide format with each IlmnID as a column and each row as a sample
  methylation_wide <- merged_data %>%
   select(IlmnID, Full_Sentrix, MethylationValue, braaksc) %>%  # Select only necessary columns
    pivot_wider(
     names_from = IlmnID,            # Each IlmnID becomes a column
     values_from = MethylationValue  # Fill values with MethylationValue
   )

  # Remove rows with missing Braak staging, if any
  methylation_wide <- methylation_wide %>%
   filter(!is.na(braaksc))  # Remove rows with missing Braak staging
#

  
# END OF STEP 1
#######################


  
  
  



  
#######################
##STEP 2: PERFORM GAMS WITH PROGRESS SAVING
#
#
  
# Specify Range
start_idx <- 1               
end_idx <- start_idx + 412998  
probe_subset <- data.reg.ROSMAP[start_idx:end_idx, ]
  
# Initialize an empty data frame to store results for this batch
output_file <- paste0("GAM_results_batch_", start_idx, "_to_", end_idx, ".csv")
  if (!file.exists(output_file)) {
    # Create file if it does not exist
    batch_results <- data.frame(Probe = rownames(probe_subset), p_value = NA, edf = NA, stringsAsFactors = FALSE)
    write.csv(batch_results, output_file, row.names = FALSE)
  } else {
    # Resume from existing file
    batch_results <- read.csv(output_file)
    cat("Resuming from saved progress file.\n")
  }
  
  # Ensure Braak staging values are in the same order as columns in data.reg.ROSMAP
  braaksc_values <- merged_data$braaksc[match(colnames(data.reg.ROSMAP), merged_data$Full_Sentrix)]
  
  # Loop through each probe in the subset
  for (probe in 1:nrow(probe_subset)) {
    # Skip already processed probes
    if (!is.na(batch_results$p_value[probe])) next
    
    # Extract methylation values for the current probe across all samples
    methylation_values <- as.numeric(probe_subset[probe, ])  # Convert to numeric vector
    
    # Fit the GAM model
    gam_model <- gam(braaksc_values ~ s(methylation_values, k = 4))
    
    # Extract p-value and effective degrees of freedom (edf)
    summary_gam <- summary(gam_model)
    batch_results$p_value[probe] <- summary_gam$s.pv[1]
    batch_results$edf[probe] <- summary_gam$edf[1]
    
    # Calculate and print the percentage completed
    percent_complete <- (probe / nrow(probe_subset)) * 100
    if (probe %% 100 == 0) {  # Print every 100 probes (optional to reduce console output)
      cat(sprintf("Progress: %.2f%% of probes modeled with GAM\n", percent_complete))
    }
    
    # Save progress every 500 probes or after every probe
    if (probe %% 500 == 0 || probe == nrow(probe_subset)) {
      write.csv(batch_results, output_file, row.names = FALSE)
      cat(sprintf("Progress saved at probe %d\n", probe))
    }
  }
  cat("GAM modeling complete for this batch. Results saved to:", output_file, "\n")
  
#
# END OF STEP 2
#######################
  
  
#######################
#
#STEP 3 PERFORM CROSS VALIDATION
#
  # Load necessary Braak stage values, ensuring they are numeric
  braaksc_values <- as.numeric(merged_data$braaksc[match(colnames(data.reg.ROSMAP), merged_data$Full_Sentrix)])
  
  # Load the GAM results file with progress saved (update the file path if necessary)
  gam_results <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/GAM_results_with_cv_partial.csv")
  
  # Identify probes that still need to be processed (i.e., where cv_error is NA)
  remaining_indices <- which(is.na(gam_results$cv_error))
  
  # Loop through only the remaining probes to perform cross-validation
  for (i in remaining_indices) {
    probe <- gam_results$Probe[i]
    
    # Check if the probe exists in the data
    if (probe %in% rownames(data.reg.ROSMAP)) {
      # Run cross-validation for the current probe
      gam_results$cv_error[i] <- cross_validate_gam(probe, data.reg.ROSMAP, braaksc_values)
    } else {
      cat(paste("Probe", probe, "not found in data.reg.ROSMAP. Skipping.\n"))
      next
    }
    
    # Print progress every 0.01% completion
    if (i %% (nrow(gam_results) / 10000) == 0) {
      percent_complete <- round((i / nrow(gam_results)) * 100, 4)
      cat(paste("Progress:", percent_complete, "% completed\n"))
    }
    
    # Save progress every 500 probes or after the last probe
    if (i %% 500 == 0 || i == remaining_indices[length(remaining_indices)]) {
      write.csv(gam_results, "GAM_results_with_cv_partial.csv", row.names = FALSE)
      cat(sprintf("Progress saved at probe %d\n", i))
    }
  }

  # Save the final results with cross-validation errors
  write.csv(gam_results, "GAM_results_with_cv.csv", row.names = FALSE)
  cat("Cross-validation complete! Results saved to GAM_results_with_cv.csv\n")
  #
  #Read GAM results into variable
  GAM_results_cv <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/GAM_results_with_cv.csv")
 
  #stats on the cv error
 mean(GAM_results_cv$cv_error) 
 sd(GAM_results_cv$cv_error)
 
 # Calculate Mean Squared Error (MSE)
 mse <- mean(GAM_results_cv$cv_error^2)
 cat("Mean Squared Error (MSE):", mse, "\n")
 
 # Calculate Mean Absolute Error (MAE)
 mae <- mean(abs(GAM_results_cv$cv_error))
 cat("Mean Absolute Error (MAE):", mae, "\n")
#
#END OF STEP 3
###########################
  
##########################
#
# STEP 4: MULTIPLE CORRECTIONS TESTING (FDR)
# Plot significant GAMS
 
 GAMResults <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/GAM_results_with_cv.csv")
 
 #Perform bonferroni correction
 # Apply Bonferroni correction to the p_value column
 GAMResults$adjusted_p_value <- p.adjust(GAMResults$p_value, method = "fdr")
 
 # Count the number of adjusted p-values below 0.05
 num_significant <- sum(GAMResults$adjusted_p_value < 0.1)
 
 # Print the result
 cat("Number of adjusted p-values below 0.1:", num_significant, "\n")
 
 # Set the adjusted p-value threshold
 adjusted_p_value_threshold <- 0.1
 
 # Create a new variable indicating whether each probe is significant based on adjusted p-value
 GAMResults$significant <- GAMResults$adjusted_p_value < adjusted_p_value_threshold
 
 # Filter to get only the significant probes
 significant_probes <- GAMResults %>% filter(significant == TRUE)
#
#END OF STEP 4
############################
 
############################
 # STEP 5: PLOT SIGNIFICANT GAMs
 #
 # Choose a subset of significant probes to plot (e.g., the first 10)
 subset_probes <- significant_probes$Probe[1:10]
 
 # Initialize an empty list to store plots
 gam_plots <- list()
 
 # Loop through each selected probe to fit the GAM and generate plots
 for (probe in subset_probes) {
   # Extract methylation values for the current probe across all samples
   methylation_values <- as.numeric(data.reg.ROSMAP[probe, ])
   
   # Ensure Braak staging values are aligned with the samples
   braaksc_values <- merged_data$braaksc[match(colnames(data.reg.ROSMAP), merged_data$Full_Sentrix)]
   
   # Fit the GAM model for the selected probe
   gam_model <- gam(methylation_values ~ s(braaksc_values, k = 4))
   
   # Create a data frame for plotting
   plot_data <- data.frame(
     Braak_Stage = braaksc_values,
     Methylation = methylation_values
   )
   
   # Generate the plot with ggplot2
   p <- ggplot(plot_data, aes(x = Methylation, y = Braak_Stage)) +
     geom_point(color = "blue", alpha = 0.5) +
     stat_smooth(method = "gam", formula = y ~ s(x, k = 4), color = "red", linewidth = 1) +
     labs(
       title = paste("GAM Fit for Probe", probe),
       x = "Methylation Level",
       y = "Braak Stage"
     ) +
     theme_minimal()
   
   # Save each plot in the list
   gam_plots[[probe]] <- p
 }
 
 # Display the plots
 for (p in gam_plots) {
   print(p)
 }
 
#
# END OF STEP 5
###############################
  
  
###############################
#
# STEP 6 : ENSEMBLE FEATURE SELECTION
 # Load necessary libraries

 
 # Step 6.1: Prepare the Data
 # Match significant probes with their methylation data
 selected_probes <- significant_probes$Probe
 rf_data_subset <- data.reg.ROSMAP[selected_probes, ]
 braaksc_values <- merged_data$braaksc[match(colnames(rf_data_subset), merged_data$Full_Sentrix)]
 
 #Normalise function
 normalize <- function(x) {
   if (max(x) == min(x)) return(rep(0, length(x)))
   (x - min(x)) / (max(x) - min(x))
 }
 
 # Normalize rf_data_subset using the normalize function
 rf_data_df <- as.data.frame(t(rf_data_subset))  # Transpose for modeling
 rf_data_df <- as.data.frame(lapply(rf_data_df, normalize)) # Ensure numeric for regression
 rf_data_normalized <- rf_data_df
 
 
 
 # Step 6.2: Train Random Forest with Cross-Validation
 train_control <- trainControl(
   method = "cv",
   number = 5,
   verboseIter = TRUE
 )
 
 # Fit Random Forest
 rf_model <- caret::train(
   x = rf_data_df,
   y = braaksc_values,
   method = "rf",
   trControl = train_control,
   tuneGrid = expand.grid(mtry = sqrt(ncol(rf_data_df))),
   ntree = 500
 )
 
 # Extract Random Forest feature importance
 rf_importance <- varImp(rf_model, scale = TRUE)$importance
 rf_importance_normalized <- normalize(rf_importance$Overall)  # Normalize RF importance
 rf_sorted_importance <- rf_importance[order(-rf_importance$Overall), , drop = FALSE]
 
 # Step 6.3: Compute Spearman Correlation Importance
 spearman_importance <- sapply(
   rf_data_df,
   function(x) cor(x, braaksc_values, method = "spearman")
 )
 spearman_importance <- abs(spearman_importance)  # Take absolute values
 spearman_importance_normalized <- normalize(spearman_importance)
 
 # Step 6.4: Fit Elastic Net Regression
 elastic_net_model <- cv.glmnet(
   as.matrix(rf_data_df),
   braaksc_values,
   alpha = 0,  # Mix between Ridge (0) and Lasso (1)
   family = "gaussian"
 )
 
 # Extract Elastic Net coefficients
 elastic_net_importance <- abs(as.vector(coef(elastic_net_model, s = "lambda.min")))[-1]  # Exclude intercept
 names(elastic_net_importance) <- colnames(rf_data_df)
 elastic_net_importance_normalized <- normalize(elastic_net_importance)
 
 # Step 6.5: Fit Gradient Boosting Machine (GBM)
 gbm_model <- gbm(
   formula = braaksc_values ~ ., 
   data = data.frame(rf_data_df, braaksc_values = braaksc_values),
   distribution = "gaussian",
   n.trees = 500,
   interaction.depth = 3,
   shrinkage = 0.01,
   cv.folds = 5
 )
 
 gbm_importance <- summary(gbm_model, plotit = FALSE)$rel.inf
 names(gbm_importance) <- colnames(rf_data_df)
 gbm_importance_normalized <- normalize(gbm_importance)
 
 
 # Step 6.6: Combine Importance Scores
 # Ensure that `rf_importance_normalized` is a named vector that matches rownames
 rf_importance_vector <- as.vector(rf_importance_normalized)
 names(rf_importance_vector) <- rownames(rf_sorted_importance)
 
 # Combine all importance scores into a single data frame
 combined_importance <- data.frame(
   Probe = rownames(rf_sorted_importance),
   RFImportance = rf_importance_vector,  # Reference the named vector
   SpearmanImportance = spearman_importance_normalized[rownames(rf_sorted_importance)],
   ElasticNetImportance = elastic_net_importance_normalized[rownames(rf_sorted_importance)],
   GBMImportance = gbm_importance_normalized[rownames(rf_sorted_importance)]
 )
 
 # Calculate the average importance, excluding `NA`
 combined_importance$AverageImportance <- rowMeans(
   combined_importance[, c("RFImportance", "SpearmanImportance", "ElasticNetImportance", "GBMImportance")],
   na.rm = TRUE
 )
 
 # Sort by Average Importance
 sorted_combined_importance <- combined_importance %>%
   arrange(desc(AverageImportance))
 
 # Save the combined importance to a CSV
 write.csv(sorted_combined_importance, "combined_feature_importance.csv", row.names = FALSE)
 
 # Visualize the top 20 features
 barplot(
   sorted_combined_importance$AverageImportance[1:20],
   names.arg = sorted_combined_importance$Probe[1:20],
   las = 2,
   main = "Top 20 Features by Average Importance"
 )
#
# END OF STEP 6
######################################
  
  
  
 ###############################
 #
 # STEP 6 : ENSEMBLE FEATURE SELECTION
 #
 ###############################
 
 # Step 6.1: Prepare the Data
 selected_probes <- significant_probes$Probe
 rf_data_subset <- data.reg.ROSMAP[selected_probes, ]
 braaksc_values <- merged_data$braaksc[match(colnames(rf_data_subset), merged_data$Full_Sentrix)]
 
 # Normalize function
 normalize <- function(x) {
   if (max(x) == min(x)) return(rep(0, length(x)))
   (x - min(x)) / (max(x) - min(x))
 }
 
 # Normalize rf_data_subset
 rf_data_df <- as.data.frame(t(rf_data_subset))  
 rf_data_df <- as.data.frame(lapply(rf_data_df, normalize)) 
 rf_data_normalized <- rf_data_df
 
 # Train Random Forest with Cross-Validation
 train_control <- trainControl(method = "cv", number = 5, verboseIter = TRUE)
 
 rf_model <- caret::train(
   x = rf_data_df,
   y = braaksc_values,
   method = "rf",
   trControl = train_control,
   tuneGrid = expand.grid(mtry = sqrt(ncol(rf_data_df))),
   ntree = 500
 )
 
 rf_importance <- varImp(rf_model, scale = TRUE)$importance
 rf_importance_normalized <- normalize(rf_importance$Overall)
 rf_sorted_importance <- rf_importance[order(-rf_importance$Overall), , drop = FALSE]
 
 # Compute Spearman Correlation Importance
 spearman_importance <- sapply(rf_data_df, function(x) cor(x, braaksc_values, method = "spearman"))
 spearman_importance <- abs(spearman_importance)  
 spearman_importance_normalized <- normalize(spearman_importance)
 
 # Fit Elastic Net Regression
 elastic_net_model <- cv.glmnet(
   as.matrix(rf_data_df),
   braaksc_values,
   alpha = 0.5,  
   family = "gaussian"
 )
 
 elastic_net_importance <- abs(as.vector(coef(elastic_net_model, s = "lambda.min")))[-1]
 names(elastic_net_importance) <- colnames(rf_data_df)
 elastic_net_importance_normalized <- normalize(elastic_net_importance)
 
 # Fit Gradient Boosting Machine (GBM)
 gbm_model <- gbm(
   formula = braaksc_values ~ ., 
   data = data.frame(rf_data_df, braaksc_values = braaksc_values),
   distribution = "gaussian",
   n.trees = 500,
   interaction.depth = 3,
   shrinkage = 0.01,
   cv.folds = 5
 )
 
 gbm_importance <- summary(gbm_model, plotit = FALSE)$rel.inf
 names(gbm_importance) <- colnames(rf_data_df)
 gbm_importance_normalized <- normalize(gbm_importance)
 
 # Combine Importance Scores
 rf_importance_vector <- as.vector(rf_importance_normalized)
 names(rf_importance_vector) <- rownames(rf_sorted_importance)
 
 combined_importance <- data.frame(
   Probe = rownames(rf_sorted_importance),
   RFImportance = rf_importance_vector,
   SpearmanImportance = spearman_importance_normalized[rownames(rf_sorted_importance)],
   ElasticNetImportance = elastic_net_importance_normalized[rownames(rf_sorted_importance)],
   GBMImportance = gbm_importance_normalized[rownames(rf_sorted_importance)]
 )
 
 combined_importance$AverageImportance <- rowMeans(
   combined_importance[, c("RFImportance", "SpearmanImportance", "ElasticNetImportance", "GBMImportance")],
   na.rm = TRUE
 )
 
 sorted_combined_importance <- combined_importance %>%
   arrange(desc(AverageImportance))
 
 write.csv(sorted_combined_importance, "combined_feature_importance.csv", row.names = FALSE)
 
 barplot(
   sorted_combined_importance$AverageImportance[1:20],
   names.arg = sorted_combined_importance$Probe[1:20],
   las = 2,
   main = "Top 20 Features by Average Importance"
 )
 
 #########################################
 #
 # STEP 7: GENE SET ENRICHMENT ANALYSIS (GSEA) using fGSEA
 #
 #########################################
 
 # Load feature selection results
 feature_importance <- read.csv("combined_feature_importance.csv")
 
 # Filter significant probes
 top_probes <- feature_importance %>%
   filter(AverageImportance > 0.5) %>%  
   select(Probe)
 
 # Merge with probe-to-gene mapping
 significant_genes <- probes %>%
   filter(IlmnID %in% top_probes$Probe) %>%
   select(UCSC_RefGene_Name) %>%
   distinct()
 
 gene_list <- as.character(significant_genes$UCSC_RefGene_Name)
 
 # Convert Gene Symbols to Entrez IDs
 entrez_ids <- bitr(
   gene_list,
   fromType = "SYMBOL",
   toType = "ENTREZID",
   OrgDb = org.Hs.eg.db
 )
 entrez_ids <- na.omit(entrez_ids)
 
 # Prepare ranked gene list for fGSEA
 ranked_genes <- feature_importance %>%
   filter(Probe %in% top_probes$Probe) %>%
   arrange(desc(AverageImportance)) %>%
   select(Probe, AverageImportance)
 
 ranked_genes$Gene <- probes$UCSC_RefGene_Name[match(ranked_genes$Probe, probes$IlmnID)]
 ranked_genes <- ranked_genes[!is.na(ranked_genes$Gene), ]
 ranks <- setNames(ranked_genes$AverageImportance, ranked_genes$Gene)
 
 # Load KEGG pathways
 kegg_pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
 
 # Run fGSEA
 fgsea_results <- fgsea(
   pathways = split(kegg_pathways$gene_symbol, kegg_pathways$gs_name),
   stats = ranks,
   minSize = 10,
   maxSize = 500,
   nperm = 1000
 )
 
 # Filter significant pathways (FDR < 0.05)
 significant_pathways <- fgsea_results %>%
   filter(padj < 0.05) %>%
   arrange(padj)
 
 # Save GSEA results
 write.csv(significant_pathways, "GSEA_results.csv", row.names = FALSE)
 
 # Visualize top pathways
 plotEnrichment(
   split(kegg_pathways$gene_symbol, kegg_pathways$gs_name)[[significant_pathways$pathway[1]]],
   ranks
 ) + labs(title = significant_pathways$pathway[1])
 
 cat("GSEA completed. Results saved in GSEA_results.csv\n")
#
# END OF STEP 7
# END OF SCRIPT
########################################################################################
  
  
  
  
  

  

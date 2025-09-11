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

source("cross_validate_gam.R")

# Set working directory (update path as needed)
setwd("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP")

# Load / Read necessary files
pheno <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/ROSMAPphenometh.csv")
samples <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/samplesROSMAP.csv")
load("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/regressedROSMAP.rda")
probes <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/450kinfo.csv")
#
#
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
##STEP 2:  PERFORM GAMS
#
#
  #Specify Range
  start_idx <- 1               
  end_idx <- start_idx + 412998  
  probe_subset <- data.reg.ROSMAP[start_idx:end_idx, ]

  # Initialize an empty data frame to store results for this batch
  batch_results <- data.frame(Probe = rownames(probe_subset), p_value = NA, edf = NA, stringsAsFactors = FALSE)

  # Ensure Braak staging values are in the same order as columns in data.reg.ROSMAP
  braaksc_values <- merged_data$braaksc[match(colnames(data.reg.ROSMAP), merged_data$Full_Sentrix)]

  # Loop through each probe in the subset
  for (probe in 1:nrow(probe_subset)) {
   # Extract methylation values for the current probe across all samples
   methylation_values <- as.numeric(probe_subset[probe, ])  # Convert to numeric vector
  
   # Fit the GAM model
   gam_model <- gam(methylation_values ~ s(braaksc_values, k = 4))
  
   # Extract p-value and effective degrees of freedom (edf)
   summary_gam <- summary(gam_model)
   batch_results$p_value[probe] <- summary_gam$s.pv[1]
   batch_results$edf[probe] <- summary_gam$edf[1]
  
   # Calculate and print the percentage completed
   percent_complete <- (probe / nrow(probe_subset)) * 100
   if (probe %% 100 == 0) {  # Print every 100 probes (optional to reduce console output)
     cat(sprintf("Progress: %.2f%% of probes modeled with GAM\n", percent_complete))
   }
  }
  # Save or print results
  final_gam_results <- batch_results
  write.csv(batch_results, paste0("GAM_results_batch_", start_idx, "_to_", end_idx, ".csv"), row.names = FALSE)
#
#
# END OF STEP 2



##########################
#STEP 3: CROSS VALIDATE GAM
#
  braaksc_values <- merged_data$braaksc[match(colnames(data.reg.ROSMAP), merged_data$Full_Sentrix)]

  # Initialize a column for cross-validation errors
  final_gam_results$cv_error <- NA
  #CALL CROSS VALIDATION FUNCTION
  for (probe in final_gam_results$Probe) {
    final_gam_results$cv_error[final_gam_results$Probe == probe] <- cross_validate_gam(probe, data.reg.ROSMAP, braaksc_values)
  }
  # Save the updated results
  write.csv(final_gam_results, "final_gam_results_with_cv.csv", row.names = FALSE)
#
#END OF STEP 3
########################


###########################
#STEP 4: PLOT SIGNIFICANT GAMs
#
#
#
#Save GAM results as variable
  final_gam_results <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/GAM_results2.csv")

# Set the p-value threshold
  p_value_threshold <- 0.05

# Create a new variable indicating whether each probe is significant
  final_gam_results$significant <- final_gam_results$p_value < p_value_threshold

# Filter to get only the significant probes
  significant_probes <- final_gam_results %>% filter(significant == TRUE)

# View the significant probes
  head(significant_probes)
#
############################


###########################
#STEP 4: PLOT SIGNIFICANT GAMs
#
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
      stat_smooth(method = "gam", formula = y ~ s(x, k = 3), color = "red", linewidth = 1) +
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
#
#END OF STEP 4
###############################

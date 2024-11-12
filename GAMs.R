# Load necessary libraries
library(dplyr)
library(tidyr)
library(randomForest)
library(mgcv)
library(ggplot2)
library(xgboost)





# Set working directory (update path as needed)
setwd("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP")

# Load / Read necessary files
pheno <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/ROSMAPphenometh.csv")
samples <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/samplesROSMAP.csv")
load("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/regressedROSMAP.rda")
probes <- read.csv("C:/Users/at883/OneDrive - University of Exeter/Desktop/ROSMAP/450kinfo.csv")

  
  #ORGANISE THE DATA
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



#PERFORM GAMS
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


#LOAD GAM RESULTS AND FILTER FOR SIGNIFICANCE
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



#PLOT
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


#RANDOM FOREST
# Extract probe names as a vector from the "Probe" column
significant_probe_names <- significant_probes$Probe




# GRADIENT BOOSTING
# Filter `merged_data` to include only significant probes
significant_probe_names <- significant_probes$Probe
gbm_data <- merged_data %>%
  filter(IlmnID %in% significant_probe_names) %>%
  select(Full_Sentrix, IlmnID, MethylationValue, braaksc)

# Convert data to a wide format with each probe as a column
gbm_wide_data <- gbm_data %>%
  pivot_wider(names_from = IlmnID, values_from = MethylationValue) %>%
  distinct()

# Separate Braak staging and predictors
gbm_labels <- as.numeric(factor(gbm_wide_data$braaksc)) - 1
gbm_matrix <- as.matrix(gbm_wide_data %>% select(-Full_Sentrix, -braaksc))

# Convert to DMatrix for XGBoost
dtrain <- xgb.DMatrix(data = gbm_matrix, label = gbm_labels)

# Set XGBoost parameters
params <- list(
  objective = "multi:softmax",
  num_class = length(unique(gbm_labels)),
  eta = 0.01,
  max_depth = 5,
  subsample = 0.8,
  lambda = 2,
  colsample_bytree = 0.8
)

# Train the Gradient Boosting model
nrounds <- 50
gbm_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = nrounds,
  watchlist = list(train = dtrain),
  verbose = 1
)

# Evaluate model performance
preds <- predict(gbm_model, dtrain)
pred_labels <- as.integer(preds)
confusion_matrix <- table(Predicted = pred_labels, Actual = gbm_labels)
print(confusion_matrix)

# Feature importance and top 5,000 probes
importance_matrix <- xgb.importance(model = gbm_model)
top_5000_probes <- importance_matrix %>%
  arrange(desc(Gain)) %>%
  head(5000) %>%
  pull(Feature)

# RANDOM FOREST with Top 5,000 Probes
# Subset `data.reg.ROSMAP` to include only the top 5,000 probes
rf_data <- data.reg.ROSMAP[rownames(data.reg.ROSMAP) %in% top_5000_probes, ]

# Filter `merged_data` to include only rows that match top probes
filtered_braaksc <- merged_data %>%
  filter(IlmnID %in% top_5000_probes) %>%
  select(Full_Sentrix, braaksc) %>%
  distinct()

# Align `filtered_braaksc` with columns in `rf_data`
braaksc_values <- filtered_braaksc$braaksc[match(colnames(rf_data), filtered_braaksc$Full_Sentrix)]
valid_indices <- !is.na(braaksc_values)
rf_data <- rf_data[, valid_indices]
braaksc_values <- factor(braaksc_values[valid_indices])

# Run the Random Forest model
rf_model <- randomForest(
  x = t(rf_data),
  y = braaksc_values,
  ntree = 50000,
  do.trace = 500,
  nodesize = 10,
  mtry = min(5000, ncol(rf_data))  # mtry is now set dynamically based on columns
)

# Print Random Forest model summary
print(rf_model)
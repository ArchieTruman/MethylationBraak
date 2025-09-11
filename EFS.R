## =========================
## STEP 6: Ensemble Feature Selection (RF, Spearman, Elastic Net, GBM)
## =========================

suppressPackageStartupMessages({
  library(caret)
  library(glmnet)
  library(gbm)
  library(dplyr)
  library(randomForest)
})

# ---- Choose dataset (EC or PFC) ----
# Example: run for EC (use LBBA_reg + braak_A for PFC)
beta_mat  <- LBBE_reg
braak_vec <- braak_E

# ---- Select significant probes from STEP 4 ----
selected_probes <- intersect(significant_probes$Probe, rownames(beta_mat))
cat("Using", length(selected_probes), "significant probes for ensemble feature selection.\n")

# Subset and transpose to samples x probes
rf_data_subset <- beta_mat[selected_probes, , drop = FALSE]
rf_data_df <- as.data.frame(t(rf_data_subset))  # rows = samples, cols = probes

# ---- Normalization function ----
normalize <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

rf_data_df <- as.data.frame(lapply(rf_data_df, normalize)) # normalize each probe

# ---- Step 6.2: Random Forest ----
set.seed(42)
train_control <- trainControl(method = "cv", number = 5, verboseIter = TRUE)

rf_model <- caret::train(
  x = rf_data_df,
  y = braak_vec,
  method = "rf",
  trControl = train_control,
  tuneGrid = expand.grid(mtry = floor(sqrt(ncol(rf_data_df)))),
  ntree = 500
)

rf_importance <- varImp(rf_model, scale = TRUE)$importance
rf_importance_vector <- normalize(rf_importance$Overall)
names(rf_importance_vector) <- rownames(rf_importance)

# ---- Step 6.3: Spearman correlation ----
spearman_importance <- sapply(rf_data_df, function(x) cor(x, braak_vec, method = "spearman"))
spearman_importance <- abs(spearman_importance)
spearman_importance_normalized <- normalize(spearman_importance)

# ---- Step 6.4: Elastic Net ----
elastic_net_model <- cv.glmnet(
  as.matrix(rf_data_df),
  braak_vec,
  alpha = 0.5,   # mix Ridge/Lasso; adjust as needed
  family = "gaussian"
)

elastic_net_importance <- abs(as.vector(coef(elastic_net_model, s = "lambda.min")))[-1] # drop intercept
names(elastic_net_importance) <- colnames(rf_data_df)
elastic_net_importance_normalized <- normalize(elastic_net_importance)

# ---- Step 6.5: Gradient Boosting ----
gbm_model <- gbm(
  formula = braak_vec ~ .,
  data = data.frame(rf_data_df, braak_vec = braak_vec),
  distribution = "gaussian",
  n.trees = 500,
  interaction.depth = 3,
  shrinkage = 0.01,
  cv.folds = 5,
  verbose = FALSE
)

gbm_importance <- summary(gbm_model, plotit = FALSE)$rel.inf
names(gbm_importance) <- colnames(rf_data_df)
gbm_importance_normalized <- normalize(gbm_importance)

# ---- Step 6.6: Combine importance scores ----
all_probes <- colnames(rf_data_df)

combined_importance <- data.frame(
  Probe = all_probes,
  RFImportance = rf_importance_vector[all_probes],
  SpearmanImportance = spearman_importance_normalized[all_probes],
  ElasticNetImportance = elastic_net_importance_normalized[all_probes],
  GBMImportance = gbm_importance_normalized[all_probes],
  stringsAsFactors = FALSE
)

combined_importance$AverageImportance <- rowMeans(
  combined_importance[, c("RFImportance", "SpearmanImportance", "ElasticNetImportance", "GBMImportance")],
  na.rm = TRUE
)

# ---- Sort and save ----
sorted_combined_importance <- combined_importance %>%
  arrange(desc(AverageImportance))

write.csv(sorted_combined_importance, "EC_combined_feature_importance.csv", row.names = FALSE)

# ---- Plot top 20 ----
barplot(
  sorted_combined_importance$AverageImportance[1:20],
  names.arg = sorted_combined_importance$Probe[1:20],
  las = 2,
  main = "Top 20 Features by Average Importance",
  cex.names = 0.7
)

cat("STEP 6 complete: results written to EC_combined_feature_importance.csv\n")



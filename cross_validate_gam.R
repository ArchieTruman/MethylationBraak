cross_validate_gam <- function(probe, methylation_data, braaksc_values, k_folds = 5) {
  # Extract methylation values for the current probe
  methylation_values <- as.numeric(methylation_data[probe, ])
  
  # Fit the GAM model
  gam_model <- gam(methylation_values ~ s(braaksc_values, k = 4))
  
  # Perform cross-validation
  cv_error <- cv.glm(data = data.frame(braaksc_values, methylation_values), glmfit = gam_model, K = k_folds)$delta[1]
  return(cv_error)
}
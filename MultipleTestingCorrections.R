## =========================
## STEP 4: Multiple testing correction (FDR + Bonferroni)
## =========================

# ---- Load GAM results ----
# Replace with whichever results you want to correct (EC or PFC, or merge them)
GAMResults <- read.csv("EC_GAM_results_1_to_230375.csv")  # example file from STEP 2
# Or: GAMResults <- read.csv("PFC_GAM_results_1_to_n.csv")

## ---- Multiple testing: Benjamini–Hochberg (FDR) ----
GAMResults$BH_p <- p.adjust(GAMResults$p_value, method = "BH")  # or method = "fdr"

# ---- Threshold ----
fdr_thresh <- 0.10

# ---- Count significant ----
num_sig_bh <- sum(GAMResults$BH_p < fdr_thresh, na.rm = TRUE)
cat("BH (FDR)-significant probes (<", fdr_thresh, "):", num_sig_bh, "\n")

# ---- Flag significance ----
GAMResults$significant <- GAMResults$BH_p < fdr_thresh

# ---- Subset significant ----
significant_probes <- GAMResults %>% dplyr::filter(significant)

# ---- Quick diagnostic plot ----
hist(GAMResults$BH_p, breaks = 50,
     main = "BH (FDR)-adjusted p-values",
     xlab = "Adjusted p-value")
abline(v = fdr_thresh, col = "red", lwd = 2)

# ---- Save outputs ----
write.csv(GAMResults, "EC_GAM_results_with_corrections.csv", row.names = FALSE)
write.csv(significant_probes, "EC_GAM_significant_probes.csv", row.names = FALSE)

cat("STEP 4 complete. BH(FDR)-corrected results written to disk.\n")

cat("STEP 4 complete. Corrected results written to disk.\n")

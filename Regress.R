LBBE <- load("DataLBBEPFdasen.rda")
PhenoLBBEPF <- read.csv("Epheno.csv")
LBBA <- load("LBBApfdasen.rda")
PhenoLBBA <- read.csv("Apheno.csv")

head(PhenoLBBEPF)
head(PhenoLBBA)

## ---------- Helpers ----------
residualize_with_covars <- function(beta_mat, pheno, formula_rhs) {
  # Build design matrix from pheno using a formula like "~ Age + Sex + prop"
  X <- model.matrix(as.formula(paste("~", formula_rhs)), data = pheno)
  
  resids <- matrix(NA_real_, nrow = nrow(beta_mat), ncol = ncol(beta_mat),
                   dimnames = dimnames(beta_mat))
  
  for (i in seq_len(nrow(beta_mat))) {
    y  <- beta_mat[i, ]
    ok <- is.finite(y) & stats::complete.cases(X)
    if (sum(ok) >= ncol(X)) {
      fit <- lm.fit(x = X[ok, , drop = FALSE], y = y[ok])
      res <- y
      res[ok] <- y[ok] - as.vector(X[ok, , drop = FALSE] %*% fit$coefficients)
      resids[i, ] <- res
    } else {
      resids[i, ] <- NA_real_
    }
  }
  resids
}

prep_beta <- function(df) {
  stopifnot(names(df)[1] == "X")     # first col is probe ID
  mat <- as.matrix(df[, -1, drop=FALSE])
  rownames(mat) <- df$X
  # strip leading "X" that R added to sample IDs
  colnames(mat) <- sub("^X", "", colnames(mat))
  mat
}

align_pheno <- function(pheno, sample_ids, id_col = "Row.names") {
  # ensure the pheno rownames are the sample IDs (no leading X)
  ids <- sub("^X", "", pheno[[id_col]])
  rownames(pheno) <- ids
  pheno[sample_ids, , drop = FALSE]
}

## ---------- Choose how to encode SEX ----------
# Option A (recommended): use the string 'Gender.x' and coerce to factor
sex_from_strings <- function(ph) factor(ph$Gender.x)

# Option B: use numeric 'Gender2' and label levels
sex_from_numeric <- function(ph) factor(ph$Gender2, levels = c(1, 2), labels = c("MALE", "FEMALE"))

## ---------- 1) EC set: LBBE + PhenoLBBEPF ----------
LBBE_beta <- prep_beta(LBBE)
PhenoLBBEPF_aligned <- align_pheno(PhenoLBBEPF, colnames(LBBE_beta), id_col = "Row.names")

# Build covariates: Age, Sex, prop
PhenoLBBEPF_aligned$Age <- as.numeric(PhenoLBBEPF_aligned$Age..at.death..x)
PhenoLBBEPF_aligned$Sex <- sex_from_strings(PhenoLBBEPF_aligned)  # or sex_from_numeric(...)
PhenoLBBEPF_aligned$prop <- as.numeric(PhenoLBBEPF_aligned$prop)

# Quick sanity checks
stopifnot(!any(is.na(PhenoLBBEPF_aligned$Age)))
stopifnot(all(levels(PhenoLBBEPF_aligned$Sex) %in% c("MALE", "FEMALE")))
stopifnot(is.numeric(PhenoLBBEPF_aligned$prop))

# Residualise: Age + Sex + prop
LBBE.reg <- residualize_with_covars(LBBE_beta, PhenoLBBEPF_aligned, formula_rhs = "Age + Sex + prop")
write.csv(LBBE.reg, "LBBE_regressed.csv", quote = FALSE)

## ---------- 2) PFC set: LBBA + PhenoLBBA ----------
LBBA_beta <- prep_beta(LBBA)
PhenoLBBA_aligned <- align_pheno(PhenoLBBA, colnames(LBBA_beta), id_col = "Row.names")

PhenoLBBA_aligned$Age  <- as.numeric(PhenoLBBA_aligned$Age..at.death..x)
PhenoLBBA_aligned$Sex  <- sex_from_strings(PhenoLBBA_aligned)  # or sex_from_numeric(...)
PhenoLBBA_aligned$prop <- as.numeric(PhenoLBBA_aligned$prop)

stopifnot(!any(is.na(PhenoLBBA_aligned$Age)))
stopifnot(all(levels(PhenoLBBA_aligned$Sex) %in% c("MALE", "FEMALE")))
stopifnot(is.numeric(PhenoLBBA_aligned$prop))

LBBA.reg <- residualize_with_covars(LBBA_beta, PhenoLBBA_aligned, formula_rhs = "Age + Sex + prop")
write.csv(LBBA.reg, "LBBA_regressed.csv", quote = FALSE)

cat("Done. Wrote LBBE_regressed.csv and LBBA_regressed.csv\n")
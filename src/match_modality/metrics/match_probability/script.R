cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/match_modality/test_resource.test_sol.h5ad",
  input_prediction = "resources_test/match_modality/test_resource.prediction.h5ad",
  output = "resources_test/match_modality/test_resource.scores.h5ad"
)
## VIASH END

cat("Read solution h5ad\n")
ad_sol <- anndata::read_h5ad(par$input_solution, backed = TRUE)

cat("Read prediction h5ad\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)

cat("Unscrambling predictions\n")
pairing_ix <- ad_sol$uns[["pairing_ix"]]
X_pred <- as(ad_pred$X, "CsparseMatrix")[, order(pairing_ix)]
dimnames(X_pred) <- NULL

cat("Calculating normalisation factors\n")
rowSum <- Matrix::rowSums(X_pred)
colSum <- Matrix::colSums(X_pred)

cat("Computing the match modality score\n")
match_probability_per_mod1 <- diag(X_pred) / rowSum
match_probability_per_mod2 <- diag(X_pred) / colSum

match_probability_mod1 <- mean(match_probability_per_mod1)
match_probability_mod2 <- mean(match_probability_per_mod2)

cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = c("match_probability_mod1", "match_probability_mod2"),
    metric_values = c(match_probability_mod1, match_probability_mod2),
    per_cell = list(
      match_probability_per_mod1 = match_probability_per_mod1, 
      match_probability_per_mod2 = match_probability_per_mod2
    )
  )
)

# should we also save the metrics object?
# this would allow for plotting the auroc and aupr curves afterwards.

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
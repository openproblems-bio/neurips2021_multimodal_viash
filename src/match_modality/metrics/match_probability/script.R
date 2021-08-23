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
ad_sol <- anndata::read_h5ad(par$input_solution)

cat("Read prediction h5ad\n")
expect_true(
  grepl("\\.h5ad$", par$input_prediction),
  info = "Prediction file should be an h5ad file"
)
ad_pred <-
  tryCatch({
    anndata::read_h5ad(par$input_prediction)
  }, error = function(e) {
    stop(paste0("Can't open prediction h5ad file. Detailed error message:\n", e$message))
  })
expect_true(
  ad_sol$uns$dataset_id == ad_pred$uns$dataset_id
)
X_sol <- ad_sol$X
X_pred <- as(ad_pred$X, "CsparseMatrix")
dimnames(X_sol) <- dimnames(X_pred) <- list(NULL, NULL)

cat("Calculating the match modality score")
X_pred_normrow <- X_pred
normrow_factor <- 1 / Matrix::rowSums(X_pred)
X_pred_normrow@x <- X_pred_normrow@x * normrow_factor[X_pred@i+1]

X_pred_normcol <- X_pred
normcol_factor <- 1 / Matrix::colSums(X_pred)
normcol_fac <- unlist(lapply(seq_len(nrow(X_pred)), function(i) rep(normcol_factor[[i]], X_pred@p[[i+1]] - X_pred@p[[i]])))
X_pred_normcol@x <- X_pred_normcol@x * normcol_fac

match_probability_mod1 <- sum(X_pred_normrow * X_sol) / sum(X_sol)
match_probability_mod2 <- sum(X_pred_normcol * X_sol) / sum(X_sol)

cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = c("match_probability_mod1", "match_probability_mod2"),
    metric_values = c(match_probability_mod1, match_probability_mod2)
  )
)

# should we also save the metrics object?
# this would allow for plotting the auroc and aupr curves afterwards.

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
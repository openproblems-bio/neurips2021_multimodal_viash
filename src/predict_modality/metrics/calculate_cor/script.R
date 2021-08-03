cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = c("resources_test/predict_modality/test_resource.solution.h5ad"),
  input_prediction = c("resources_test/predict_modality/test_resource.prediction.h5ad"),
  output = "test_resource.scores.h5ad"
)
# par <- list(
#   input_solution = list.files("/tmp/neurips2021_work/51/422eb9e5c4c90c15083880ef26a94f", full.names = TRUE, recursive = TRUE, pattern = "*.output_solution.h5ad"),
#   input_prediction = list.files("/tmp/neurips2021_work/51/422eb9e5c4c90c15083880ef26a94f", full.names = TRUE, recursive = TRUE, pattern = "*.output.h5ad"),
#   output = "test_resource.scores.h5ad"
# )
## VIASH END

cat("Reading solution file\n")
ad_sol <- 
  tryCatch({
    anndata::read_h5ad(par$input_solution)
  }, error = function(e) {
    cat("Warning! Could not read solution file '", sol_path, "'. Error:\n", sep = "")
    cat(e$message)
    NULL
  })

cat("Reading prediction file\n")
ad_pred <- 
  tryCatch({
    anndata::read_h5ad(par$input_prediction)
  }, error = function(e) {
    cat("Warning! Could not read prediction file '", sol_path, "'. Error:\n", sep = "")
    cat(e$message)
    NULL
  })

cat("Check prediction format\n")
expect_equal(
  ad_sol$uns$dataset_id, ad_pred$uns$dataset_id,
  info = "Prediction and solution have differing dataset_ids"
)
expect_true(
  all.equal(dim(ad_sol), dim(ad_pred)),
  info = "Dataset and prediction anndata objects should have the same shape / dimensions."
)

cat("Computing correlation metrics\n")
# Wrangle data
tv <- ad_sol$X
pv <- ad_pred$X

# Compute metrics
pearson_mat <- dynutils::calculate_similarity(tv, pv, method = "pearson", margin = 2, diag = TRUE, drop0 = TRUE)
spearman_mat <- dynutils::calculate_similarity(tv, pv, method = "spearman", margin = 2, diag = TRUE, drop0 = TRUE)
score_pearson <- mean(diag(pearson_mat)) / 2 + 0.5
score_spearman <- mean(diag(spearman_mat)) / 2 + 0.5

cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = c("score_pearson", "score_spearman"),
    metric_values = c(score_pearson, score_spearman),
    metric_moreisbetter = c(TRUE, TRUE)
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

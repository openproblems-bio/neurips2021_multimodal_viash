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
#   input_solution = list.files("/tmp/neurips2021_work/62/df32515759063107a687584c24a303", full.names = TRUE, recursive = TRUE, pattern = "*.output_solution.h5ad"),
#   input_prediction = list.files("/tmp/neurips2021_work/62/df32515759063107a687584c24a303", full.names = TRUE, recursive = TRUE, pattern = "*.output.h5ad"),
#   output = "test_resource.scores.h5ad"
# )
## VIASH END

cat("Reading solution file\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

cat("Reading prediction file\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)

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

mean_pearson <- mean(diag(pearson_mat))
mean_spearman <- mean(diag(spearman_mat))
score_mean_pearson <- mean_pearson / 2 + .5
score_mean_spearman <- mean_spearman / 2 + .5

metric_ids <- c("mean_pearson", "mean_spearman", "score_mean_pearson", "score_mean_spearman")
metric_values <- c(mean_pearson, mean_spearman, score_mean_pearson, score_mean_spearman)

cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = metric_ids,
    metric_values = metric_values
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

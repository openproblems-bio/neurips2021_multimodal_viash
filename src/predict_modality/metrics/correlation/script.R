cat("Load dependencies\n")
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/predict_modality/test_resource.test_mod2.h5ad",
  input_prediction = "resources_test/predict_modality/test_resource.prediction.h5ad",
  output = "test_resource.scores.h5ad"
)
par <- list(
  input_solution = "work/40/ef654fcb79dac281fc9d856eb773e6/totalvi_10x_pbmc_10k_rna.censor_dataset.output_test_mod2.h5ad",
  input_prediction = "work/40/ef654fcb79dac281fc9d856eb773e6/totalvi_10x_pbmc_10k_rna.dummy_random.output.h5ad",
  output = "test_resource.scores.h5ad"
)
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
pearson_vec <- diag(dynutils::calculate_similarity(tv, pv, method = "pearson", margin = 2, diag = TRUE, drop0 = TRUE))
spearman_vec <- diag(dynutils::calculate_similarity(tv, pv, method = "spearman", margin = 2, diag = TRUE, drop0 = TRUE))

pearson_vec[!is.finite(pearson_vec)] <- 0
spearman_vec[!is.finite(spearman_vec)] <- 0

mean_pearson <- mean(pearson_vec)
mean_spearman <- mean(spearman_vec)

metric_ids <- c("mean_pearson", "mean_spearman")
metric_values <- c(mean_pearson, mean_spearman)

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

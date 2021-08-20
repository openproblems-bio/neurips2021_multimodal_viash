cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = c("resources_test/predict_modality/test_resource.test_mod2.h5ad"),
  input_prediction = c("resources_test/predict_modality/test_resource.prediction.h5ad"),
  output = "test_resource.scores.h5ad"
)
## VIASH END

cat("Reading solution file\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

cat("Reading prediction file\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)
# ad_pred$X <- as(ad_pred$X * 0 + 1, "CsparseMatrix")

cat("Check prediction format\n")
expect_equal(
  ad_sol$uns$dataset_id, ad_pred$uns$dataset_id,
  info = "Prediction and solution have differing dataset_ids"
)
expect_true(
  all.equal(dim(ad_sol), dim(ad_pred)),
  info = "Dataset and prediction anndata objects should have the same shape / dimensions."
)

cat("Computing MSE and MSLE metrics\n")
# Wrangle data
scores <- 
  full_join(
    summary(ad_sol$X) %>% rename(solx = x),
    summary(ad_pred$X) %>% rename(predx = x),
    by = c("i", "j")
  ) %>%
  select(-i, -j) %>%
  mutate(
    solx = ifelse(is.na(solx), 0, solx),
    predx = ifelse(is.na(predx), 0, predx),
  ) %>%
  summarise(
    mse = mean((solx - predx)^2),
    msle = mean((log(solx + 1) - log(abs(predx) + 1) * sign(predx))^2)
  )

cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = names(scores),
    metric_values = as.vector(as.matrix(scores))
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

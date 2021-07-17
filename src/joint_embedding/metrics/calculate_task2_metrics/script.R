cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/task2/test_resource.solution.h5ad",
  input_prediction = "resources_test/task2/test_resource.prediction.h5ad",
  output = "resources_test/task2/test_resource.scores.h5ad"
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

cat("Calculating metrics\n")
df <- data.frame(ad_pred$X, SOLUTION_CELL_TYPE = ad_sol$obs[["cell_type"]])
rf <- ranger::ranger(SOLUTION_CELL_TYPE ~ ., df)

rf_oob_correct_pred <- 1 - rf$prediction.error


cat("Create output object\n")
out <- anndata::AnnData(
  X = NULL,
  shape = dim(ad_sol),
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = c("rf_oob_correct_pred"),
    metric_values = c(rf_oob_correct_pred),
    metric_moreisbetter = c(TRUE)
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

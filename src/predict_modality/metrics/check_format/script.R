cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
task <- "predict_modality"
par <- list(
  input_solution = paste0("resources_test/", task, "/test_resource.test_mod2.h5ad"),
  input_prediction = paste0("resources_test/", task, "/test_resource.prediction.h5ad"),
  output = paste0("resources_test/", task, "/test_resource.scores.h5ad")
)
## VIASH END

cat("Read prediction h5ad\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

correct_api <- tryCatch({
  # read prediction
  ad_pred <- anndata::read_h5ad(par$input_prediction)

  # check dataset id
  dataset_id <- ad_pred$uns[["dataset_id"]]
  assert_that(dataset_id == ad_sol$uns[["dataset_id"]])

  # check method id
  method_id <- ad_pred$uns[["method_id"]]
  assert_that(
    is.character(method_id),
    method_id != ""
  )

  # check X
  X <- ad_pred$X
  expected_X <- ad_sol$X
  assert_that(
    is(X, "sparseMatrix"),
    nrow(X) == nrow(expected_X),
    ncol(X) == ncol(expected_X),
    !is.null(rownames(X)),
    !is.null(colnames(X)),
    all(rownames(X) == rownames(expected_X)),
    all(colnames(X) == colnames(expected_X))
  )

  1
}, error = function(e) {
  cat("ERROR: ", e$message, "\n", sep = "")
  0
})


cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = c("finished", "correct_format"),
    metric_values = c(1, correct_format)
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

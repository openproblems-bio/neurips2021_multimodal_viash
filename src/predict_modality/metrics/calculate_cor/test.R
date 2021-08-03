library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

par <- list(
  input_solution = c("resources_test/predict_modality/test_resource.solution.h5ad"),
  input_prediction = c("resources_test/predict_modality/test_resource.prediction.h5ad"),
  output = "test_resource.scores.h5ad"
)

cat("> Running metrics\n")
out <- processx::run(
  command = "./calculate_cor",
  args = c(
    "--input_solution", par$input_solution,
    "--input_prediction", par$input_prediction,
    "--output", par$output
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists(par$output))

cat("> Checking contents of output.h5ad\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)
ad_out <- anndata::read_h5ad(par$output)

expect_equal(ad_out$uns[["dataset_id"]], ad_pred$uns[["dataset_id"]])
expect_equal(ad_out$uns[["method_id"]], ad_pred$uns[["method_id"]])
expect_gte(length(ad_out$uns[["metric_ids"]]), 1)
expect_equal(length(ad_out$uns[["metric_ids"]]), length(ad_out$uns[["metric_values"]]))

cat("> Test succeeded!\n")

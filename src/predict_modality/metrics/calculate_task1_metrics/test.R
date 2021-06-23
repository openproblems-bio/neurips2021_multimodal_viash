library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

assert_that(
  file.exists("dataset_task1_solution.h5ad"),
  file.exists("dataset_task1_prediction_randomforest.h5ad")
)

cat("> Running metrics\n")
out <- processx::run(
  command = "./calculate_task1_metrics",
  args = c(
    "--input_solution", "dataset_task1_solution.h5ad",
    "--input_prediction", "dataset_task1_prediction_randomforest.h5ad",
    "--output", "dataset_task1_prediction_randomforest_scores.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(
  file.exists("dataset_task1_prediction_randomforest_scores.h5ad")
)

cat("> Checking contents of output.h5ad\n")
adata_orig <- anndata::read_h5ad("dataset_task1_prediction_randomforest.h5ad")
adata <- anndata::read_h5ad("dataset_task1_prediction_randomforest_scores.h5ad")

assert_that(
  adata$uns[["dataset_id"]] == adata_orig$uns[["dataset_id"]],
  adata$uns[["method_id"]] == adata_orig$uns[["method_id"]],
  adata$n_obs == adata_orig$n_obs,
  adata$n_vars == adata_orig$n_vars,
  length(adata$uns[["metric_ids"]]) == length(adata$uns[["metric_values"]])
)

# TODO: check content of layers["prediction"]

cat("> Test succeeded!\n")

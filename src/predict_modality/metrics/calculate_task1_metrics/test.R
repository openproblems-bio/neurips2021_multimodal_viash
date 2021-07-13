library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

par <- list(
  input_solution = "resources_test/task1/pbmc_1k_protein_v3.solution.h5ad",
  input_prediction = "resources_test/task1/pbmc_1k_protein_v3.prediction.h5ad",
  output = "resources_test/task1/pbmc_1k_protein_v3.scores.h5ad"
)

cat("> Running metrics\n")
out <- processx::run(
  command = "./calculate_task1_metrics",
  args = c(
    "--input_solution", "resources_test/task1/pbmc_1k_protein_v3.solution.h5ad",
    "--input_prediction", "resources_test/task1/pbmc_1k_protein_v3.prediction.h5ad",
    "--output", "output.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists("output.h5ad"))

cat("> Checking contents of output.h5ad\n")
adata_orig <- anndata::read_h5ad("resources_test/task1/pbmc_1k_protein_v3.prediction.h5ad")
adata <- anndata::read_h5ad("output.h5ad")

expect_equal(adata$uns[["dataset_id"]], adata_orig$uns[["dataset_id"]])
expect_equal(adata$uns[["method_id"]], adata_orig$uns[["method_id"]])
expect_gte(length(adata$uns[["metric_ids"]]), 1)
expect_equal(length(adata$uns[["metric_ids"]]), length(adata$uns[["metric_values"]]))

cat("> Test succeeded!\n")

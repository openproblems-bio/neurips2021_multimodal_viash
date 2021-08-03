library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

par <- list(
  input_solution = list.files("resources_test/joint_embedding", full.names = TRUE, recursive = TRUE, pattern = "*.solution.h5ad"),
  input_prediction = list.files("resources_test/joint_embedding", full.names = TRUE, recursive = TRUE, pattern = "*.prediction.h5ad"),
  output = "pairing.tsv"
)

# todo: create incorrect input files

cat("> Running metrics\n")
out <- processx::run(
  command = "./gatekeep_predictions",
  args = c(
    "--input_solution", paste(par$input_solution, collapse = ":"),
    "--input_prediction", paste(par$input_prediction, collapse = ":"),
    "--output", par$output
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists(par$output))

cat("> Checking contents of pairing.tsv\n")
pairing <- readr::read_tsv(par$output)

# expect_equal(adata$uns[["dataset_id"]], adata_orig$uns[["dataset_id"]])
# expect_equal(adata$uns[["method_id"]], adata_orig$uns[["method_id"]])
# expect_gte(length(adata$uns[["metric_ids"]]), 1)
# expect_equal(length(adata$uns[["metric_ids"]]), length(adata$uns[["metric_values"]]))

cat("> Test succeeded!\n")

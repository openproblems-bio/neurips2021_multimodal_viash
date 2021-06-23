library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running baseline component\n")
out <- processx::run(
  command = "./baseline_randomforest",
  args = c(
    "--input", "dataset_task1_censored.h5ad",
    "--output", "output.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(
  file.exists("output.h5ad")
)

cat("> Checking contents of output.h5ad\n")
adata_orig <- anndata::read_h5ad("dataset_task1_censored.h5ad")
adata <- anndata::read_h5ad("output.h5ad")

assert_that(
  adata$uns[["dataset_id"]] == adata_orig$uns[["dataset_id"]],
  adata$uns[["method_id"]] == "baseline_randomforest",
  adata$n_obs == adata_orig$n_obs,
  adata$n_vars == adata_orig$n_vars,
  all(c("prediction") %in% names(adata$layers))
)

# TODO: check content of layers["prediction"]

cat("> Test succeeded!\n")

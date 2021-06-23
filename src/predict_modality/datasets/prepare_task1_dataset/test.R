library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running censor component\n")
out <- processx::run(
  command = "./prepare_task1_dataset",
  args = c(
    "--input", "dataset.h5ad",
    "--output_censored", "output_censored.h5ad",
    "--output_solution", "output_solution.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(
  file.exists("output_censored.h5ad")
)

cat("> Checking contents of output_censored.h5ad\n")
adata_orig <- anndata::read_h5ad("dataset.h5ad")
adata_censor <- anndata::read_h5ad("output_censored.h5ad")

assert_that(
  adata_censor$uns[["dataset_id"]] == paste0(adata_orig$uns[["dataset_id"]], "_task1"),
  adata_censor$n_obs == adata_orig$n_obs,
  adata_censor$n_vars == adata_orig$n_vars,
  all(c("modality1", "modality2") %in% names(adata_censor$layers))
)

cat("> Checking contents of output_solution.h5ad\n")
adata_sol <- anndata::read_h5ad("output_solution.h5ad")
assert_that(
  adata_sol$uns[["dataset_id"]] == paste0(adata_orig$uns[["dataset_id"]], "_task1"),
  adata_sol$n_obs == adata_orig$n_obs,
  adata_sol$n_vars == adata_orig$n_vars,
  all(c("modality2") %in% names(adata_sol$layers))
)

# TODO: check content of layers["modality1"] and layers["modality2"]

cat("> Test succeeded!\n")

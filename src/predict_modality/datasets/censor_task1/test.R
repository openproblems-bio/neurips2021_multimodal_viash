library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running censor component\n")
out <- processx::run(
  command = "./censor_task1",
  args = c(
    "--input", "dataset.h5ad",
    "--output", "output.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(
  file.exists("output.h5ad")
)

cat("> Checking contents of output.h5ad\n")
adata_orig <- anndata::read_h5ad("dataset.h5ad")
adata <- anndata::read_h5ad("output.h5ad")

assert_that(
  adata$uns[["dataset_id"]] == paste0(adata_orig$uns[["dataset_id"]], "_task1"),
  adata$n_obs == adata_orig$n_obs,
  adata$n_vars == adata_orig$n_vars,
  all(c("modality1", "modality2") %in% names(adata$layers))
)

cat("> Test succeeded!\n")

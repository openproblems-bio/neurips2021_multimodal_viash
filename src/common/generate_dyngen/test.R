library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running dyngen\n")
out <- processx::run(
  # command = "viash",
  command = "./generate_dyngen",
  args = c(
    # "run", "src/common/generate_dyngen/config.vsh.yaml", "--",
    "--id", "mytest",
    "--backbone", "bifurcating",
    "--output", "dataset.h5ad",
    "--output_censored", "dataset_censored.h5ad",
    "--plot", "plot.pdf",
    "--model", "model.rds",
    "--num_threads", "1",
    "--num_simulations", "3",
    "--num_cells", "100",
    "--num_genes", "50",
    "--census_interval", "2",
    "--compute_rna_velocity",
    "--compute_cellwise_grn"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(file.exists("dataset.h5ad"))
assert_that(file.exists("dataset_censored.h5ad"))
assert_that(file.exists("plot.pdf"))
assert_that(file.exists("model.rds"))

cat("> Checking contents of dataset.h5ad\n")
adata <- anndata::read_h5ad("dataset.h5ad")

assert_that(
  adata$uns[["dataset_id"]] == "mytest",
  adata$n_obs == 100,
  adata$n_vars == 50,
  "regulatory_network_sc" %in% adata$obsm_keys(),
  all(c("counts_spliced", "counts_unspliced", "counts_protein") %in% names(adata$layers))
)


cat("> Checking contents of dataset_censored.h5ad\n")
adata2 <- anndata::read_h5ad("dataset_censored.h5ad")

assert_that(
  adata2$uns[["dataset_id"]] == "mytest",
  adata2$n_obs == 100,
  adata2$n_vars == 50,
  all(c("modality1", "modality2") %in% names(adata2$layers)),
  all(names(adata2$layers) %in% c("modality1", "modality2"))
)

cat("> Test succeeded!\n")

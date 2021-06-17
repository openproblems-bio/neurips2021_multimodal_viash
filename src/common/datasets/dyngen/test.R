library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running dyngen\n")
out <- processx::run(
  # command = "viash",
  command = "./dyngen",
  args = c(
    # "run", "src/common/generate_dyngen/config.vsh.yaml", "--",
    "--id", "mytest",
    "--backbone", "bifurcating",
    "--output", "dataset.h5ad",
    "--plot", "plot.pdf",
    "--num_threads", "1",
    "--num_simulations", "3",
    "--num_cells", "100",
    "--num_genes", "50",
    "--census_interval", "2",
    "--store_rna_velocity",
    "--store_atac",
    "--store_antibody"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(
  file.exists("dataset.h5ad"),
  file.exists("plot.pdf")
)

cat("> Checking contents of dataset.h5ad\n")
adata <- anndata::read_h5ad("dataset.h5ad")

assert_that(
  adata$uns[["dataset_id"]] == "mytest",
  adata$n_obs == 100,
  adata$n_vars == 50,
  all(c("antibody", "atac") %in% names(adata$layers))
)

cat("> Test succeeded!\n")

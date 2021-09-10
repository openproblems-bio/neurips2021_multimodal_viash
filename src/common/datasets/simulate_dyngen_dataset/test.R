library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running dyngen\n")
out <- processx::run(
  command = "./simulate_dyngen_dataset",
  args = c(
    "--id", "mytest",
    "--output_rna", "dataset_rna.h5ad",
    "--output_mod2", "dataset_mod2.h5ad",
    "--plot", "plot.pdf",
    "--num_threads", "1",
    "--num_simulations", "3",
    "--num_cells", "100",
    "--num_genes", "50",
    "--census_interval", "2",
    "--store_protein",
    "--num_proteins", "20"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists("dataset_rna.h5ad"))
expect_true(file.exists("dataset_mod2.h5ad"))
expect_true(file.exists("plot.pdf"))

cat("> Reading output files\n")
ad1 <- anndata::read_h5ad("dataset_rna.h5ad")
ad2 <- anndata::read_h5ad("dataset_mod2.h5ad")

cat("> Checking contents of dataset_rna.h5ad\n")
expect_equal(ad1$uns[["dataset_id"]], "mytest")
expect_equal(ad1$n_obs, 100)
expect_equal(ad1$n_vars, 50)

cat("> Checking contents of dataset_mod2.h5ad\n")
expect_equal(ad2$uns[["dataset_id"]], "mytest")
expect_equal(ad2$n_obs, 100)
expect_equal(ad2$n_vars, 20)

cat("> Test succeeded!\n")

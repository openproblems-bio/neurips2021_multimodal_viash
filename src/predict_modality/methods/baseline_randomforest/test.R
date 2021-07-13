library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running baseline component\n")
out <- processx::run(
  command = "./baseline_randomforest",
  args = c(
    "--input_mod1", "resources_test/task1/pbmc_1k_protein_v3.mod1.h5ad",
    "--input_mod2", "resources_test/task1/pbmc_1k_protein_v3.mod2.h5ad",
    "--output", "output.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists("output.h5ad"))

cat("> Checking contents of output.h5ad\n")
adata_mod1 <- anndata::read_h5ad("resources_test/task1/pbmc_1k_protein_v3.mod1.h5ad")
adata_mod2 <- anndata::read_h5ad("resources_test/task1/pbmc_1k_protein_v3.mod2.h5ad")
adata <- anndata::read_h5ad("output.h5ad")

expect_equal(adata$uns[["dataset_id"]], adata_mod1$uns[["dataset_id"]])
expect_equal(adata$uns[["method_id"]], "baseline_randomforest")
expect_equal(adata$n_obs, adata_mod1$n_obs - adata_mod2$n_obs)
expect_equal(adata$n_vars, adata_mod2$n_vars)

cat("> Test succeeded!\n")

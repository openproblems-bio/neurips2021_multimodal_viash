library(testthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running censor component\n")
out <- processx::run(
  command = "./prepare_task1_dataset",
  args = c(
    "--input_rna", "resources_test/common/pbmc_1k_protein_v3.normalize.output_rna.h5ad",
    "--input_mod2", "resources_test/common/pbmc_1k_protein_v3.normalize.output_mod2.h5ad",
    "--output_mod1", "output_mod1.h5ad",
    "--output_mod2", "output_mod2.h5ad",
    "--output_solution", "output_solution.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists("output_mod1.h5ad"))
expect_true(file.exists("output_mod2.h5ad"))
expect_true(file.exists("output_solution.h5ad"))

cat("> Reading h5ad files\n")
adata_orig <- anndata::read_h5ad("resources_test/common/pbmc_1k_protein_v3.normalize.output_rna.h5ad")
adata_mod1 <- anndata::read_h5ad("output_mod1.h5ad")
adata_mod2 <- anndata::read_h5ad("output_mod2.h5ad")
adata_sol <- anndata::read_h5ad("output_solution.h5ad")

cat("> Checking contents of h5ad files\n")
expect_equal(adata_mod1$uns[["dataset_id"]], paste0(adata_orig$uns[["dataset_id"]], "_task1"))
expect_equal(adata_mod2$uns[["dataset_id"]], paste0(adata_orig$uns[["dataset_id"]], "_task1"))
expect_equal(adata_sol$uns[["dataset_id"]], paste0(adata_orig$uns[["dataset_id"]], "_task1"))
expect_equal(adata_mod1$n_obs, adata_orig$n_obs)
expect_equal(adata_mod2$n_obs + adata_sol$n_obs, adata_orig$n_obs)

# TODO check contents of matrices, check rownames

cat("> Test succeeded!\n")

library(testthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

cat("> Running censor component\n")
out <- processx::run(
  command = "./censor_dataset",
  args = c(
    "--input_mod1", "resources_test/common/test_resource.output_rna.h5ad",
    "--input_mod2", "resources_test/common/test_resource.output_mod2.h5ad",
    "--output_train_mod1", "output_train_mod1.h5ad",
    "--output_train_mod2", "output_train_mod2.h5ad",
    "--output_test_mod1", "output_test_mod1.h5ad",
    "--output_test_mod2", "output_test_mod2.h5ad"
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists("output_train_mod1.h5ad"))
expect_true(file.exists("output_train_mod2.h5ad"))
expect_true(file.exists("output_test_mod1.h5ad"))
expect_true(file.exists("output_test_mod2.h5ad"))

cat("> Reading h5ad files\n")
adata_orig <- anndata::read_h5ad("resources_test/common/test_resource.output_rna.h5ad")
adata_train_mod1 <- anndata::read_h5ad("output_train_mod1.h5ad")
adata_train_mod2 <- anndata::read_h5ad("output_train_mod2.h5ad")
adata_test_mod1 <- anndata::read_h5ad("output_test_mod1.h5ad")
adata_test_mod2 <- anndata::read_h5ad("output_test_mod2.h5ad")

cat("> Checking contents of h5ad files\n")
expected_did <- paste0(adata_orig$uns[["dataset_id"]], "_PM_gex2adt")
expect_equal(adata_train_mod1$uns[["dataset_id"]], expected_did)
expect_equal(adata_train_mod2$uns[["dataset_id"]], expected_did)
expect_equal(adata_test_mod1$uns[["dataset_id"]], expected_did)
expect_equal(adata_test_mod2$uns[["dataset_id"]], expected_did)
expect_equal(adata_train_mod1$n_obs + adata_test_mod1$n_obs, adata_orig$n_obs)
expect_equal(adata_train_mod2$n_obs + adata_test_mod2$n_obs, adata_orig$n_obs)

# TODO check contents of matrices, check rownames

cat("> Test succeeded!\n")

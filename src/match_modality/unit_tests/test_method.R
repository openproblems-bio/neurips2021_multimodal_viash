library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
# This code block will be replaced by viash at runtime.
meta <- list(functionality_name = "foo")
## VIASH END

method_id <- meta[["functionality_name"]]
command <- paste0("./", method_id)

# define some filenames
testpar <- list(
  input_train_mod1 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad",
  input_train_mod2 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad",
  input_train_sol = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_sol.h5ad",
  input_test_mod1 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad",
  input_test_mod2 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
  input_test_sol = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_sol.h5ad",
  output = "output.h5ad"
)

cat("> Running method\n")
out <- processx::run(
  command = paste0("./", meta[["functionality_name"]]),
  args = c(
    "--input_train_mod1", testpar$input_train_mod1,
    "--input_train_mod2", testpar$input_train_mod2,
    "--input_train_sol", testpar$input_train_sol,
    "--input_test_mod1", testpar$input_test_mod1,
    "--input_test_mod2", testpar$input_test_mod2,
    "--output", testpar$output
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(file.exists(testpar$output))

cat("> Reading h5ad files\n")
ad_sol <- anndata::read_h5ad(testpar$input_test_sol)
ad_pred <- anndata::read_h5ad(testpar$output)

cat("> Checking dataset id\n")
dataset_id <- ad_pred$uns[["dataset_id"]]
assert_that(dataset_id == ad_sol$uns[["dataset_id"]])

cat("> Checking method id\n")
assert_that(ad_pred$uns[["method_id"]] == method_id)

cat("> Checking X\n")
assert_that(
  is(ad_pred$X, "sparseMatrix"),
  length(ad_pred$X@x) <= 1000 * ad_sol$n_obs,
  ad_pred$n_obs == ad_sol$n_obs,
  ad_pred$n_vars == ad_sol$n_vars,
  all(ad_pred$X@x >= 0),
  isTRUE(all.equal(
    Matrix::rowSums(ad_pred$X),
    rep(1, ad_pred$n_obs),
    check.attributes = FALSE,
    tolerance = 1e-5
  ))
)

cat("> Test succeeded!\n")

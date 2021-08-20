library(assertthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

par <- list(
  input_train_mod1 = "resources_test/match_modality/test_resource.train_mod1.h5ad",
  input_train_mod2 = "resources_test/match_modality/test_resource.train_mod2.h5ad",
  input_train_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  input_test_mod1 = "resources_test/match_modality/test_resource.test_mod1.h5ad",
  input_test_mod2 = "resources_test/match_modality/test_resource.test_mod2.h5ad",
  input_test_sol = "resources_test/match_modality/test_resource.test_sol.h5ad",
  output = "output.h5ad"
)

cat("> Running method\n")
out <- processx::run(
  command = paste0("./", meta[["functionality_name"]]),
  args = c(
    "--input_train_mod1", par$input_train_mod1,
    "--input_train_mod2", par$input_train_mod2,
    "--input_train_sol", par$input_train_sol,
    "--input_test_mod1", par$input_test_mod1,
    "--input_test_mod2", par$input_test_mod2,
    "--output", par$output
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
assert_that(file.exists(par$output))

cat("> Reading h5ad files\n")
ad_sol <- anndata::read_h5ad(par$input_test_sol)
ad_pred <- anndata::read_h5ad(par$output)

cat("> Checking dataset id\n")
dataset_id <- ad_pred$uns[["dataset_id"]]
assert_that(dataset_id == ad_sol$uns[["dataset_id"]])

cat("> Checking method id\n")
method_id <- ad_pred$uns[["method_id"]]
assert_that(
  is.character(method_id),
  method_id == meta[["functionality_name"]]
)

cat("> Checking X\n")
assert_that(
  is(ad_pred$X, "sparseMatrix"),
  ad_pred$n_obs == ad_sol$n_obs,
  ad_pred$n_vars == ad_sol$n_vars,
  !is.null(ad_pred$obs_names),
  all(ad_pred$obs_names == ad_sol$obs_names),
  !is.null(ad_pred$var_names),
  all(ad_pred$var_names == ad_sol$var_names)
)

cat("> Test succeeded!\n")

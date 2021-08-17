library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

par <- list(
  input_mod1 = "resources_test/joint_embedding/test_resource.mod1.h5ad",
  input_mod2 = "resources_test/joint_embedding/test_resource.mod2.h5ad",
  input_solution = "resources_test/joint_embedding/test_resource.solution.h5ad",
  output = "output.h5ad"
)

cat("> Running method\n")
out <- processx::run(
  command = paste0("./", meta[["functionality_name"]]),
  args = c(
    "--input_solution", par$input_solution,
    "--output", par$output
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists(par$output))

cat("> Reading h5ad files\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
output <- anndata::read_h5ad(par$output)

cat("> Checking contents of output.h5ad\n")
expect_equal(output$uns[["dataset_id"]], input_mod1$uns[["dataset_id"]])
expect_equal(output$uns[["method_id"]], meta[["functionality_name"]])
expect_equal(output$n_obs, input_mod1$n_obs)
expect_gte(output$n_vars, 1)
expect_lte(output$n_vars, 100)
expect_false(is.null(output$obs_names))
expect_equal(output$obs_names, input_mod1$obs_names)

cat("> Test succeeded!\n")

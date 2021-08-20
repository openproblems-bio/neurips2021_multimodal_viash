library(testthat, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
# This code block will be replaced by viash at runtime.
meta <- list(functionality_name = "foo")
## VIASH END

method_id <- meta[["functionality_name"]]
command <- paste0("./", method_id)

# define some filenames
testpar <- list(
  input_mod1 = "resources_test/joint_embedding/test_resource.mod1.h5ad",
  input_mod2 = "resources_test/joint_embedding/test_resource.mod2.h5ad",
  output = "output.h5ad"
)

cat("> Running method\n")
out <- processx::run(
  command = paste0("./", meta[["functionality_name"]]),
  args = c(
    "--input_mod1", testpar$input_mod1,
    "--input_mod2", testpar$input_mod2,
    "--output", testpar$output
  ),
  stderr_to_stdout = TRUE
)

cat("> Checking whether output files were created\n")
expect_true(file.exists(testpar$output))

cat("> Reading h5ad files\n")
input_mod1 <- anndata::read_h5ad(testpar$input_mod1)
output <- anndata::read_h5ad(testpar$output)

cat("> Checking contents of output.h5ad\n")
expect_equal(output$uns[["dataset_id"]], input_mod1$uns[["dataset_id"]])
expect_equal(output$uns[["method_id"]], method_id)
expect_equal(output$n_obs, input_mod1$n_obs)
expect_gte(output$n_vars, 1)
expect_lte(output$n_vars, 100)
expect_false(is.null(output$obs_names))
expect_equal(output$obs_names, input_mod1$obs_names)

cat("> Test succeeded!\n")

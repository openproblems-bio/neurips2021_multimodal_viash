cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_test_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  output = "output.h5ad",
)
## VIASH END

cat("Reading h5ad files\n")
input_test_sol <- anndata::read_h5ad(par$input_test_sol)

input_test_sol$uns[["method_id"]] <- "dummy_solution"

cat("Writing predictions to file\n")
zzz <- input_test_sol$write_h5ad(par$output, compression = "gzip")

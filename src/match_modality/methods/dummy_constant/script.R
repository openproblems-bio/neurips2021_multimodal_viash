cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_train_mod1 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad",
  input_train_mod2 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad",
  input_train_sol = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_sol.h5ad",
  input_test_mod1 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad",
  input_test_mod2 = "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
  output = "output.h5ad"
)
meta <- list(functionality_name = "foo")
## VIASH END

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
# input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
# input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
# input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1, backed = TRUE)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2, backed = TRUE)

knn_df <-
  expand.grid(
    i = seq_len(nrow(input_test_mod1)),
    j = seq_len(min(nrow(input_test_mod2), 1000))
  )
knn_mat <- 
  Matrix::sparseMatrix(
    i = knn_df$i,
    j = knn_df$j,
    x = rep(1, nrow(knn_df)),
    dims = list(nrow(input_test_mod1), nrow(input_test_mod2))
  )

# normalise to make rows sum to 1
rs <- Matrix::rowSums(knn_mat)
knn_mat@x <- knn_mat@x / rs[knn_mat@i + 1]

cat("Creating output anndata\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = input_test_mod1$uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

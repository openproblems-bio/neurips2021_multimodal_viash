cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(keras, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  # input_train_mod1 = "resources_test/match_modality/test_resource.train_mod1.h5ad",
  # input_train_mod2 = "resources_test/match_modality/test_resource.train_mod2.h5ad",
  # input_train_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  # input_test_mod1 = "resources_test/match_modality/test_resource.test_mod1.h5ad",
  # input_test_mod2 = "resources_test/match_modality/test_resource.test_mod2.h5ad",
  input_test_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  output = "output.h5ad",
)
## VIASH END

cat("Reading h5ad files\n")
# input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
# input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
# input_train_sol <- anndata::read_h5ad(par$input_train_sol)
# input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
# input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)
input_test_sol <- anndata::read_h5ad(par$input_test_sol)

cat("Writing predictions to file\n")
zzz <- input_test_sol$write_h5ad(par$output, compression = "gzip")

cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/task2/test_resource.mod1.h5ad",
  input_mod2 = "resources_test/task2/test_resource.mod2.h5ad",
  output = "output.h5ad",
  n_dims = 10L,
  n_neighbors = 15L,
  metric = "euclidean",
  n_pcs = 50L
)
## VIASH END

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

x <- cbind(ad1$X, ad2$X)

x_pca <- irlba::prcomp_irlba(
  cbind(ad1$X, ad2$X),
  n = par$n_pcs
)$x

cat("Performing DR\n")
dr <- uwot::umap(
  x_pca,
  n_components = par$n_dims,
  n_neighbors = par$n_neighbors,
  metric = par$metric,
  n_threads = 1,
  nn_method = "annoy"
)

rownames(dr) <- rownames(ad1)
colnames(dr) <- paste0("comp_", seq_len(par$n_dims))

out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = "baseline_pca"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

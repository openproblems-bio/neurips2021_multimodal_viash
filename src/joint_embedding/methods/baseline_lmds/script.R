cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
# path <- "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
path <- "output/datasets/joint_embedding/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.censor_dataset.output_"
par <- list(
  input_mod1 = paste0(path, "mod1.h5ad"),
  input_mod2 = paste0(path, "mod2.h5ad"),
  output = "output.h5ad",
  n_dims = 10L,
  distance_method = "spearman"
)
## VIASH END

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

cat("Performing DR\n")
dr <- lmds::lmds(
  cbind(ad1$X, ad2$X),
  ndim = par$n_dims,
  distance_method = par$distance_method
)

rownames(dr) <- rownames(ad1)
colnames(dr) <- paste0("comp_", seq_len(par$n_dims))

out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = "baseline_lmds"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

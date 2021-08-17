# Dependencies:
#   python: anndata
#   r: anndata, lmds
#
# R starter kit for the NeurIPS 2021 Single-Cell Competition. Parts
# with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/r/

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(lmds, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
# Anything within this block will be removed by viash
# and will be replaced with the parameters as specified in
# your config.vsh.yaml.
par <- list(
  input_mod1 = "sample_data/test_resource.mod1.h5ad",
  input_mod2 = "sample_data/test_resource.mod2.h5ad",
  output = "output.h5ad",
  distance_method = "spearman",
  n_pcs = 4L
)
## VIASH END

method_id <- "r_starter_kit" # fill in the name of your method here

cat("Reading h5ad files\n")
ad1 <- read_h5ad(par$input_mod1)
ad2 <- read_h5ad(par$input_mod2)

ad_merge <- anndata::concat(list(ad1, ad2), axis = 1L)

# TODO: implement own method

cat("Performing dimensionality reduction on the mod1 values\n")
# LMDS is more efficient than regular MDS because
# it does not compure a square distance matrix.
dr <- lmds(
  as(ad_merge$X, "CsparseMatrix"),
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

cat("Creating output AnnData\n")
out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")
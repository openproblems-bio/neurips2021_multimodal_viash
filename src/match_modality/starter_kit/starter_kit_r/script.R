# Dependencies:
#   python: anndata
#   r: anndata, lmds, FNN
#
# R starter kit for the NeurIPS 2021 Single-Cell Competition.
# Parts with `TODO` are supposed to be changed by you.
#
# More documentation:
#
# https://viash.io/docs/creating_components/r/

cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("lmds", quietly = TRUE)
requireNamespace("FNN", quietly = TRUE)

## VIASH START
# Anything within this block will be removed by viash
# and will be replaced with the parameters as specified in
# your config.vsh.yaml.

dataset_path <- "sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# dataset_path <- "output/datasets/match_modality/openproblems_bmmc_multiome_phase1_rna/openproblems_bmmc_multiome_phase1_rna.censor_dataset.output_"

par <- list(
  input_train_mod1 = paste0(dataset_path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(dataset_path, "train_mod2.h5ad"),
  input_train_sol = paste0(dataset_path, "train_sol.h5ad"),
  input_test_mod1 = paste0(dataset_path, "test_mod1.h5ad"),
  input_test_mod2 = paste0(dataset_path, "test_mod2.h5ad"),
  output = "output.h5ad",
  n_neighbors = 5L
)
## VIASH END

method_id <- "r_starter_kit" # fill in the name of your method here

# TODO: implement own method

# This starter kit is split up into several steps.
# * compute dimensionality reduction on [train_mod1, test_mod1] data
# * compute dimensionality reduction on [train_mod2, test_mod2] data
# * predict test_dr2 from [train_dr1, train_dr2, test_dr1]
# * calculate k nearest neighbors between test_dr2 and predicted pred_dr2
# * transform k nearest neighbors into a pairing matrix

cat("read train solution\n")
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
match_train <- input_train_sol$uns$pairing_ix + 1
rm(input_train_sol)
gc()

cat("compute dimensionality reduction on [train_mod1, test_mod1] data\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
dr_mod1 <- lmds::lmds(
  rbind(input_train_mod1$X, input_test_mod1$X),
  ndim = 10,
  distance_method = "pearson"
)
train_mod1_uns <- input_train_mod1$uns

# clear memory
rm(input_train_mod1, input_test_mod1)
gc()

cat("compute dimensionality reduction on [train_mod2, test_mod2] data\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

dr_mod2 <- lmds::lmds(
  rbind(input_train_mod2$X[order(match_train), , drop = FALSE], input_test_mod2$X),
  ndim = 10,
  distance_method = "pearson"
)
# clear memory
rm(input_train_mod2, input_test_mod2)
gc()

# split DR matrices
dr_mod1_train <- dr_mod1[seq_along(match_train), , drop = FALSE]
dr_mod2_train <- dr_mod2[seq_along(match_train), , drop = FALSE]
dr_mod1_test <- dr_mod1[-seq_along(match_train), , drop = FALSE]
dr_mod2_test <- dr_mod2[-seq_along(match_train), , drop = FALSE]


cat("predict test_dr2 from [train_dr1, train_dr2, test_dr1]\n")
dr_mod2_test_pred <- apply(dr_mod2_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_mod1_train,
    test = dr_mod1_test,
    y = yi,
    k = min(15, nrow(dr_mod1_test))
  )$pred
})

cat("calculate k nearest neighbors between test_dr2 and predicted pred_dr2\n")
knn_out <- FNN::get.knnx(
  dr_mod2_test_pred,
  dr_mod2_test,
  k = min(par$n_neighbors, nrow(dr_mod1_test))
)

cat("transform k nearest neighbors into a pairing matrix\n")
df <- tibble(
  i = as.vector(row(knn_out$nn.index)),
  j = as.vector(knn_out$nn.index),
  x = max(knn_out$nn.dist) * 2 - as.vector(knn_out$nn.dist)
)
knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  x = df$x,
  dims = list(nrow(dr_mod1_test), nrow(dr_mod2_test))
)

# normalise to make rows sum to 1
rs <- Matrix::rowSums(knn_mat)
knn_mat@x <- knn_mat@x / rs[knn_mat@i + 1]

cat("creating output anndata\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = train_mod1_uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

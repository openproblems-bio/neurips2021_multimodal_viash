# Dependencies:
#   python: anndata
#   r: anndata, lmds
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
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(keras, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
# Anything within this block will be removed by viash
# and will be replaced with the parameters as specified in
# your config.vsh.yaml.
par <- list(
  input_train_mod1 = "resources_test/match_modality/test_resource.train_mod1.h5ad",
  input_train_mod2 = "resources_test/match_modality/test_resource.train_mod2.h5ad",
  input_train_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  input_test_mod1 = "resources_test/match_modality/test_resource.test_mod1.h5ad",
  input_test_mod2 = "resources_test/match_modality/test_resource.test_mod2.h5ad",
  output = "output.h5ad",
  distance_method = "pearson"
)
## VIASH END

method_id <- "r_starter_kit" # fill in the name of your method here

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

match_train <- apply(input_train_sol$X, 1, function(x) which(x > 0)) %>% unname

# TODO: implement own method

# This starter kit is split up into several steps.
# * compute dimensionality reduction on [train_mod1, test_mod1] data
# * train regression model to predict the train_mod2 data from the dr_mod1 values
# * predict test_mod2 matrix from model and test_mod1
# * calculate k nearest neighbors between test_mod2 and predicted test_mod2
# * transform k nearest neighbors into a pairing matrix

cat("compute dimensionality reduction on [train_mod1, test_mod1] data\n")
# merge input matrices
mod1_X <- rbind(input_train_mod1$X, input_test_mod1$X)
mod2_X <- rbind(input_train_mod2$X[match_train, , drop = FALSE], input_test_mod2$X)

# perform DR
dr_x1 <- lmds::lmds(mod1_X, ndim = 10, distance_method = par$distance_method)
dr_x2 <- lmds::lmds(mod2_X, ndim = 3, distance_method = par$distance_method)

# split input matrices
dr_x1_train <- dr_x1[seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x2_train <- dr_x2[seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x1_test <- dr_x1[-seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x2_test <- dr_x2[-seq_len(nrow(input_train_mod1)), , drop = FALSE]

cat("train regression model to predict the train_mod2 data from the dr_mod1 values\n")
model <-
  keras_model_sequential() %>%
  layer_dense(100, "relu", input_shape = ncol(dr_x1)) %>%
  layer_dense(32, "relu") %>%
  layer_dense(ncol(dr_x2), "linear")

model %>% compile(
  loss = "mse",
  optimizer = "adam"
)
model %>% fit(dr_x1_train, dr_x2_train, epochs = 200, verbose = FALSE)

cat("predict test_mod2 matrix from model and test_mod1\n")
preds <- predict(model, dr_x1_test)
colnames(preds) <- colnames(dr_x2_test)


cat("calculate k nearest neighbors between test_mod2 and predicted test_mod2\n")
par_frac <- 1
knn_out <- FNN::get.knnx(
  preds,
  dr_x2_test,
  k = min(100, ceiling(par_frac * nrow(preds)))
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
  dims = list(nrow(input_test_mod1), nrow(input_test_mod2))
)

cat("write prediction output\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = input_train_mod1$uns[["dataset_id"]],
    method_id = method_id
  )
)
zzz <- out$write_h5ad(par$output, compression = "gzip")

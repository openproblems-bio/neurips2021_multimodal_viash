cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
path <- "output/public_datasets/predict_modality/openproblems_bmmc_multiome_mod2/openproblems_bmmc_multiome_mod2.censor_dataset.output_"
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  output = "output.h5ad",
  n_pcs = 4L,
  n_neighbors = 3,
  distance_method = "pearson"
)
## VIASH END

cat("Reading mod1 files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

cat("Performing DR on the mod1 values\n")
# LMDS is more efficient than regular MDS because
# it does not compure a square distance matrix.
dr <- lmds::lmds(
  rbind(input_train_mod1$X, input_test_mod1$X),
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

rm(input_train_mod1)
rm(input_test_mod1)

cat("Reading mod2 files\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)

ix <- seq_len(nrow(input_train_mod2))
dr_train <- dr[ix, , drop = FALSE]
dr_test <- dr[-ix, , drop = FALSE]
responses_train <- input_train_mod2$X

cat("Predicting for each column in modality 2\n")
preds <- apply(responses_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_train,
    test = dr_test,
    y = yi,
    k = par$n_neighbors
  )$pred
})

cat("Creating outputs object\n")
# store prediction as a sparse matrix
prediction <- Matrix::Matrix(
  preds,
  sparse = TRUE,
  dimnames = list(rownames(dr_test), colnames(input_train_mod2))
)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = input_train_mod2$uns[["dataset_id"]],
    method_id = "baseline_knn"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

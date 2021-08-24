cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_train_mod1 = "resources_test/predict_modality/test_resource.train_mod1.h5ad",
  input_test_mod1 = "resources_test/predict_modality/test_resource.test_mod1.h5ad",
  input_train_mod2 = "resources_test/predict_modality/test_resource.train_mod2.h5ad",
  output = "output.h5ad",
  n_pcs = 4L,
  n_neighbors = 3
)
## VIASH END

cat("Reading h5ad files\n")
ad1_train <- anndata::read_h5ad(par$input_train_mod1)
ad1_test <- anndata::read_h5ad(par$input_test_mod1)
ad2 <- anndata::read_h5ad(par$input_train_mod2)

cat("Performing DR on the mod1 values\n")
# LMDS is more efficient than regular MDS because
# it does not compure a square distance matrix.
dr <- lmds::lmds(
  rbind(ad1_train$X, ad1_test$X),
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

ix <- seq_len(nrow(ad1_train))
dr_train <- dr[ix, , drop = FALSE]
dr_test <- dr[-ix, , drop = FALSE]
responses_train <- ad2$X

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
  dimnames = list(rownames(dr_test), colnames(ad2))
)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = ad1_train$uns[["dataset_id"]],
    method_id = "baseline_knn"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

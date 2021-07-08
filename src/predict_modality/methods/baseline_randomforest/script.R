cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("randomForest", quietly = TRUE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/task1/pbmc_1k_protein_v3.mod1.h5ad",
  input_mod2 = "resources_test/task1/pbmc_1k_protein_v3.mod2.h5ad",
  output = "resources_test/task1/pbmc_1k_protein_v3.prediction.h5ad",
  n_pcs = 4L
)
## VIASH END

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

cat("Performing DR on the mod1 values\n")
dr <- prcomp(ad1$X, rank. = 4)$x

dr_train <- dr[ad1$obs$split == "train",]
responses_train <- ad2$X
dr_test <- dr[ad1$obs$split == "test",]

cat("Predicting for each column in modality 2\n")
preds <- lapply(seq_len(ncol(responses_train)), function(i) {
  y <- responses_train[,i]
  if (length(unique(y)) > 1) {
    rf <- randomForest::randomForest(x = dr_train, y = y)
    stats::predict(rf, dr_test)
  } else {
    setNames(rep(unique(y), nrow(dr_test)), rownames(dr_test))
  }
})

cat("Creating outputs object\n")
prediction <- Matrix::Matrix(do.call(cbind, preds), sparse = TRUE)
rownames(prediction) <- rownames(dr_test)
colnames(prediction) <- colnames(ad2)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = "baseline_randomforest"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

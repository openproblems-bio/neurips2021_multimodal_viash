cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/task1/pbmc_1k_protein_v3.mod1.h5ad",
  input_mod2 = "resources_test/task1/pbmc_1k_protein_v3.mod2.h5ad",
  output = "output.h5ad",
  n_pcs = 4L
)
## VIASH END

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

cat("Performing DR on the mod1 values\n")
dr <- lmds::lmds(
  ad1$X, 
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

dr_train <- dr[ad1$obs$split == "train",]
responses_train <- ad2$X
dr_test <- dr[ad1$obs$split == "test",]

cat("Predicting for each column in modality 2\n")
preds <- lapply(seq_len(ncol(responses_train)), function(i) {
  y <- responses_train[,i]
  yout <- 
    if (length(unique(y)) > 1) {
      FNN::knn.reg(
        train = dr_train, 
        test = dr_test,
        y = y,
        k = par$n_neighbors
      )$pred
    } else {
      rep(unique(y), nrow(dr_test))
    }
  setNames(yout, rownames(dr_test))
})

cat("Creating outputs object\n")
prediction <- Matrix::Matrix(do.call(cbind, preds), sparse = TRUE)
rownames(prediction) <- rownames(dr_test)
colnames(prediction) <- colnames(ad2)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = "baseline_knearestneighbors"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

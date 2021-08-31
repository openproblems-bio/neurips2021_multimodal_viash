cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("ranger", quietly = TRUE)

## VIASH START
par <- list(
  input_train_mod1 = "resources_test/predict_modality/test_resource.train_mod1.h5ad",
  input_test_mod1 = "resources_test/predict_modality/test_resource.test_mod1.h5ad",
  input_train_mod2 = "resources_test/predict_modality/test_resource.train_mod2.h5ad",
  output = "output.h5ad",
  n_pcs = 4L,
  distance_method = "pearson"
)
## VIASH END

cat("Reading mod1 files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

cat("Performing DR on the mod1 values\n")
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
responses_train <- input_train_mod2$X
dr_test <- dr[-ix, , drop = FALSE]

cat("Predicting for each column in modality 2\n")
preds <- lapply(seq_len(ncol(responses_train)), function(i) {
  y <- responses_train[, i]
  pred <-
   if (length(unique(y)) > 1) {
      rf <- ranger::ranger(
        data = data.frame(YTRAIN = y, dr_train),
        dependent.variable.name = "YTRAIN",
        num.threads = NULL # set to null for auto parallellisation
      )
      stats::predict(rf, dr_test)$predictions
    } else {
      rep(unique(y), nrow(dr_test))
    }
  pred
})

cat("Creating outputs object\n")
prediction <- Matrix::Matrix(do.call(cbind, preds), sparse = TRUE)
rownames(prediction) <- rownames(dr_test)
colnames(prediction) <- colnames(input_train_mod2)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = input_train_mod2$uns[["dataset_id"]],
    method_id = "baseline_randomforest"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

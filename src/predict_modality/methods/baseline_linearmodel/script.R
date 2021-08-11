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
  n_pcs = 4L
)
## VIASH END

cat("Reading h5ad files\n")
ad1_train <- anndata::read_h5ad(par$input_train_mod1)
ad1_test <- anndata::read_h5ad(par$input_test_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

cat("Performing DR on the mod1 values\n")
dr <- lmds::lmds(
  rbind(ad1_train$X, ad1_test$X), 
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

dr_train <- dr[1:nrow(ad1_train),]
responses_train <- ad2$X
dr_test <- dr[(nrow(ad1_train)+1):nrow(dr),]

cat("Predicting for each column in modality 2\n")
preds <- lapply(seq_len(ncol(responses_train)), function(i) {
  y <- responses_train[,i]
  if (length(unique(y)) > 1) {
    lm <- lm(YTRAIN~., data.frame(dr_train, YTRAIN=y))
    stats::predict(lm, data.frame(dr_test))
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
    dataset_id = ad2$uns[["dataset_id"]],
    method_id = "baseline_linearmodel"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

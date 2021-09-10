cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
path <- "output/datasets/predict_modality/openproblems_bmmc_multiome_phase1_mod2/openproblems_bmmc_multiome_phase1_mod2.censor_dataset.output_"
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

n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

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
ix <- seq_len(nrow(input_train_mod1))
dr_train <- dr[ix, , drop = FALSE]
dr_test <- dr[-ix, , drop = FALSE]

# free memory
rm(input_train_mod1, input_test_mod1)
rm()
gc()

cat("Reading mod2 files\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)

cat("Predicting for each column in modality 2\n")
pred_df <- bind_rows(pbapply::pblapply(
  seq_len(ncol(input_train_mod2)),
  cl = n_cores,
  function(j) {
    out <- FNN::knn.reg(
      train = dr_train,
      test = dr_test,
      y = input_train_mod2$X[,j],
      k = par$n_neighbors
    )$pred
    ix <- which(out != 0)
    if (length(ix) > 0) {
      tibble(
        i = ix,
        j = j,
        x = out[ix]
      )
    } else {
      NULL
    }
  }
))

cat("Creating outputs object\n")
# store prediction as a sparse matrix
prediction <- Matrix::sparseMatrix(
  i = pred_df$i,
  j = pred_df$j,
  x = pred_df$x,
  dim = c(nrow(dr_test), ncol(input_train_mod2)),
  dimnames = list(rownames(dr_test), colnames(input_train_mod2))
)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = input_train_mod2$uns[["dataset_id"]],
    method_id = meta$functionality_name
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

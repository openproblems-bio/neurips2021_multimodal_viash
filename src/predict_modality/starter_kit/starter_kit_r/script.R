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
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("lmds", quietly = TRUE)
requireNamespace("FNN", quietly = TRUE)

## VIASH START
# Anything within this block will be removed by viash
# and will be replaced with the parameters as specified in
# your config.vsh.yaml.

dataset_path <- "sample_data/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# dataset_path <- "output/datasets/predict_modality/openproblems_bmmc_multiome_phase1_rna/openproblems_bmmc_multiome_phase1_rna.censor_dataset.output_"

par <- list(
  input_train_mod1 = paste0(dataset_path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(dataset_path, "train_mod2.h5ad"),
  input_test_mod1 = paste0(dataset_path, "test_mod1.h5ad"),
  output = "output.h5ad",
  n_pcs = 4L,
  distance_method = "pearson",
  n_neighbors = 3
)
## VIASH END

method_id <- "r_starter_kit" # fill in the name of your method here

cat("Reading mod1 h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
train_mod1_uns <- input_train_mod1$uns

cat("Performing DR on the mod1 values\n")
# LMDS is more efficient than regular MDS because
# it does not compure a square distance matrix.
dr_mod1 <- lmds::lmds(
  rbind(input_train_mod1$X, input_test_mod1$X),
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

ix <- seq_len(nrow(input_train_mod1))
dr_mod1_train <- dr_mod1[ix, , drop = FALSE]
dr_mod1_test <- dr_mod1[-ix, , drop = FALSE]

# remove previous objects to save memory
rm(input_train_mod1, input_test_mod1)
gc()

cat("Reading mod2 h5ad files\n")
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)

cat("Predicting for each column in modality 2\n")
preds <- apply(input_train_mod2$X, 2, function(yi) {
  FNN::knn.reg(
    train = dr_mod1_train,
    test = dr_mod1_test,
    y = yi,
    k = par$n_neighbors
  )$pred
})

cat("Creating outputs object\n")
# store prediction as a sparse matrix
prediction <- Matrix::Matrix(
  preds,
  sparse = TRUE,
  dimnames = list(rownames(dr_mod1_test), colnames(input_train_mod2))
)

out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = train_mod1_uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

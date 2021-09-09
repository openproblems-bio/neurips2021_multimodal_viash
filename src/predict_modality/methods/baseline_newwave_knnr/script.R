cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("NewWave", quietly = TRUE)
requireNamespace("FNN", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
# path <- "resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
path <- "output/public_datasets/predict_modality/dyngen_citeseq_3_rna/dyngen_citeseq_3_rna.censor_dataset.output_"
# path <- "debug/debug."
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  output = "output.h5ad",
  newwave_maxiter = 10L,
  newwave_ngene = 100L,
  newwave_ncell = 100L
)
meta <- list(functionality_name = "foo")
## VIASH END

n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)

# fetch batch labels
batch1 <- c(as.character(input_train_mod1$obs$batch), as.character(input_test_mod1$obs$batch))
batch2 <- as.character(input_train_mod1$obs$batch)

cat("Running NewWave on mod1\n")
data1 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = cbind(t(input_train_mod1$layers[["counts"]]), t(input_test_mod1$layers[["counts"]]))),
  colData = data.frame(batch = factor(batch1))
)
data1 <- data2[Matrix::rowSums(SummarizedExperiment::assay(data1)) > 0, ]
res1 <- NewWave::newWave(
  data1,
  X = "~batch",
  verbose = TRUE,
  K = 10,
  maxiter_optimize = par$newwave_maxiter,
  n_gene_par = min(par$newwave_ngene, nrow(data1)),
  n_cell_par = min(par$newwave_ncell, ncol(data1)),
  commondispersion = FALSE
)
dr_x1 <- SingleCellExperiment::reducedDim(res1)
colnames(dr_x1) <- paste0("comp_", seq_len(ncol(dr_x1)))

# split DR matrices
train_ix <- seq_len(nrow(input_train_mod1))
dr_x1_train <- dr_x1[train_ix, , drop = FALSE]
dr_x1_test <- dr_x1[-train_ix, , drop = FALSE]

cat("Predicting mod1 DR of test cells\n")
pred <- map_df(seq_len(ncol(input_train_mod2)), function(j) {
  out <- FNN::knn.reg(
    train = dr_x1_train,
    test = dr_x1_test,
    y = input_train_mod2$X[,j],
    k = min(15, nrow(dr_x1_test))
  )$pred
  tibble(
    i = which(out > 0), 
    j = j, 
    x = out[i]
  )
})

cat("Creating outputs object\n")
# store prediction as a sparse matrix
prediction <- Matrix::sparseMatrix(
  i = pred$i,
  j = pred$j,
  x = pred$x,
  dimnames = list(rownames(input_test_mod1), colnames(input_train_mod2))
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

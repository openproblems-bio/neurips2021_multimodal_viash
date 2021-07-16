cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(lmds, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
# Anything within this block will be removed by viash
# and will be replaced with the parameters as specified in
# your config.vsh.yaml.
par <- list(
  input_mod1 = "sample_data/pbmc_1k_protein_v3.mod1.h5ad",
  input_mod2 = "sample_data/pbmc_1k_protein_v3.mod2.h5ad",
  distance_method = "spearman",
  output = "output.h5ad",
  n_pcs = 4L
)
## VIASH END

method_id = "mymethod" # fill in the name of your method here

cat("Reading h5ad files\n")
ad1 <- read_h5ad(par$input_mod1)
ad2 <- read_h5ad(par$input_mod2)

cat("Performing dimensionality reduction on the mod1 values\n")
dr <- lmds(
  ad1$X,
  ndim = par$n_pcs,
  distance_method = par$distance_method
)

# split up the train vs. test dimensionality reduction
dr_train <- dr[ad1$obs$group == "train",]
dr_test <- dr[ad1$obs$group == "test",]
responses_train <- ad2$X[,i]

cat("Run KNN regression.\n")
# For every column in mod2, predict mod2 for the test cells using the K nearest mod2 neighbors
preds <- apply(responses_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_train, 
    test = dr_test,
    y = yi,
    k = par$n_neighbors
  )$pred
})

cat("Creating output matrix\n")
prediction <- Matrix::Matrix(
  preds, 
  sparse = TRUE,
  dimnames = list(rownames(dr_test), colnames(ad2))
)

cat("Creating output AnnData\n")
out <- anndata::AnnData(
  X = prediction,
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
out$write_h5ad(par$output, compression = "gzip")

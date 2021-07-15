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

cat("For every column in mod2, train a linear model and generate predictions.\n")
preds <- lapply(seq_len(ncol(ad2)), function(i) {
  train_data <- data.frame(dr_train, predictorcolumn = ad2$X[,i])
  test_data <- data.frame(dr_test)

  # train model on train cells
  lm <- lm(predictorcolumn ~ ., train_data)

  # generate predictions on test cells
  predict(lm, test_data)
})

cat("Creating output matrix\n")
prediction <- Matrix(do.call(cbind, preds), sparse = TRUE)
rownames(prediction) <- rownames(dr_test)
colnames(prediction) <- colnames(ad2)

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

## VIASH START
par <- list(
  input = "resources/test/dyngen_bifurcating_antibody/dataset.h5ad",
  # input = "resources/test/dyngen_bifurcating_atac/dataset.h5ad",
  output = "output.h5ad"
)
## VIASH END

# load libraries
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)

# load h5ad file
adata <- anndata::read_h5ad(par$input)

# check for modalities
has_antibody <- "antibody" %in% names(adata$layers)
has_atac <- "atac" %in% names(adata$layers)

if (has_antibody == has_atac) {
  stop("Strictly one of adata.layers[\"atac\"] and adata.layers[\"antibody\"] must be defined.")
}

if (has_antibody) {
  modality1 <- adata$X
  modality2 <- adata$layers[["antibody"]]
}

if (has_atac) {
  modality1 <- adata$layers[["atac"]]
  modality2 <- adata$X
}

## TODO: should I also bother censoring the cell names and the gene names?

# detect features the second modality has
ix <- which(colSums(modality2) > 0)

# split features into (train, test)
split <- sample.int(length(ix), 0.2 * length(ix), replace = FALSE)
train_ix <- ix[-split]
test_ix <- ix[split]

# throw away 'test' features
modality2_train <- modality2
modality2_train[, test_ix] <- 0
modality2_train <- Matrix::drop0(modality2_train)

# create censored dataset
out <- anndata::AnnData(
  X = NULL,
  shape = dim(modality1),
  layers = list(
    modality1 = modality1,
    modality2 = modality2_train
  ),
  var = data.frame(
    row.names = colnames(modality2),
    is_predictor = seq_len(ncol(modality2)) %in% train_ix,
    is_response = seq_len(ncol(modality2)) %in% test_ix
  ),
  uns = list(
    dataset_id = paste0(adata$uns[["dataset_id"]], "_task1")
  )
)

# save as h5ad
zzz <- out$write_h5ad(par$output, compression = "gzip")

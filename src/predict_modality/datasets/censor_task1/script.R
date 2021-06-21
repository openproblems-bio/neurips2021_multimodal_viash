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
requireNamespace("assertthat", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)

# load h5ad file
adata <- anndata::read_h5ad(par$input)

# check for modalities
has_protein <- "protein" %in% names(adata$layers)
has_chromatin <- "chromatin" %in% names(adata$layers)

if (has_protein == has_chromatin) {
  stop("Strictly one of adata.layers[\"chromatin\"] and adata.layers[\"protein\"] must be defined.")
}

if (has_protein) {
  modality1 <- adata$X
  modality2 <- adata$layers[["protein"]]
}

if (has_chromatin) {
  modality1 <- adata$layers[["chromatin"]]
  modality2 <- adata$X
}

## TODO: should I also bother censoring the cell names and the gene names?

# throw away 'test' features
modality2_notest <- modality2
modality2_notest[adata$obs$experiment == "test", ] <- 0
modality2_notest <- Matrix::drop0(modality2_notest)

# create censored dataset
out <- anndata::AnnData(
  X = NULL,
  shape = dim(modality1),
  layers = list(
    modality1 = modality1,
    modality2 = modality2_notest
  ),
  obs = adata$obs %>% select(experiment),
  uns = list(
    dataset_id = paste0(adata$uns[["dataset_id"]], "_task1")
  )
)
assertthat::assert_that(sum(out$layers["modality2"][adata$obs$experiment == "test", ]) == 0)

# save as h5ad
zzz <- out$write_h5ad(par$output, compression = "gzip")

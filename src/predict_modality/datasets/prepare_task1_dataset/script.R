## VIASH START
par <- list(
  input = "resources/test/dyngen_bifurcating_antibody/dataset.h5ad",
  output_censored = "resources/test/dyngen_bifurcating_antibody/dataset_task1_censor.h5ad",
  output_solution = "resources/test/dyngen_bifurcating_antibody/dataset_task1_solution.h5ad"
  # input = "resources/test/dyngen_bifurcating_atac/dataset.h5ad",
)
## VIASH END

###############################################################################
###                            LOAD DEPENDENCIES                            ###
###############################################################################
# load libraries
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)


###############################################################################
###                             READ INPUT DATA                             ###
###############################################################################
# load h5ad file
adata <- anndata::read_h5ad(par$input)

# check for modalities
has_protein <- "protein" %in% names(adata$layers)
has_chromatin <- "chromatin" %in% names(adata$layers)
assert_that(
  has_protein != has_chromatin,
  msg = "Strictly one of adata.layers[\"chromatin\"] and adata.layers[\"protein\"] must be defined."
)

if (has_protein) {
  modality1 <- adata$X
  modality2 <- adata$layers[["protein"]]
  modality_types <- c(modality1 = "mrna", modality2 = "protein")
} else if (has_chromatin) {
  modality1 <- adata$layers[["chromatin"]]
  modality2 <- adata$X
  modality_types <- c(modality1 = "chromatin", modality2 = "mrna")
}


###############################################################################
###                          CREATE CENSOR OBJECT                           ###
###############################################################################

## TODO: should I also bother censoring the cell names and the gene names?

# throw away 'test' features
modality2_notest <- modality2
modality2_notest[adata$obs$experiment == "test", ] <- 0
modality2_notest <- Matrix::drop0(modality2_notest)

# create censored dataset
out_censor <- anndata::AnnData(
  X = NULL,
  shape = dim(modality1),
  layers = list(
    modality1 = modality1,
    modality2 = modality2_notest
  ),
  obs = adata$obs %>% select(experiment),
  uns = list(
    dataset_id = paste0(adata$uns[["dataset_id"]], "_task1"),
    modality_types = modality_types
  )
)
assert_that(
  sum(abs(out_censor$layers["modality2"][adata$obs$experiment == "test", ])) == 0,
  msg = "modality2 values for test cells should sum to 0."
)

# save as h5ad
zzz <- out_censor$write_h5ad(par$output_censored, compression = "gzip")


###############################################################################
###                          CREATE SOLUTION OBJECT                         ###
###############################################################################

# throw away 'test' features
modality2_onlytest <- modality2
modality2_onlytest[adata$obs$experiment != "test", ] <- 0
modality2_onlytest <- Matrix::drop0(modality2_notest)

# create censored dataset
out_solution <- anndata::AnnData(
  X = NULL,
  shape = dim(modality1),
  layers = list(
    modality2 = modality2_onlytest
  ),
  obs = adata$obs %>% select(experiment),
  uns = list(
    dataset_id = paste0(adata$uns[["dataset_id"]], "_task1"),
    modality_types = modality_types
  )
)
assert_that(
  sum(abs(out_solution$layers["modality2"][adata$obs$experiment == "test", ])) == 0,
  msg = "modality2 values for test cells should sum to 0."
)

# save as h5ad
zzz <- out_solution$write_h5ad(par$output_solution, compression = "gzip")


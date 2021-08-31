cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
par <- list(
  # input_mod1 = "resources_test/common/test_resource.output_rna.h5ad",
  # input_mod2 = "resources_test/common/test_resource.output_mod2.h5ad",
  # output_train_mod1 = "output_train_mod1.h5ad",
  # output_train_mod2 = "output_train_mod2.h5ad",
  # output_test_mod1 = "output_test_mod1.h5ad",
  # output_test_mod2 = "output_test_mod2.h5ad",
  input_mod1 = "output/public_datasets/common/openproblems_bmmc_multiome/openproblems_bmmc_multiome.manual_formatting.output_mod2.h5ad",
  input_mod2 = "output/public_datasets/common/openproblems_bmmc_multiome/openproblems_bmmc_multiome.manual_formatting.output_rna.h5ad",
  output_train_mod1 = "output/public_datasets/predict_modality/openproblems_bmmc_multiome_mod2/openproblems_bmmc_multiome_mod2.censor_dataset.output_train_mod1.h5ad",
  output_train_mod2 = "output/public_datasets/predict_modality/openproblems_bmmc_multiome_mod2/openproblems_bmmc_multiome_mod2.censor_dataset.output_train_mod2.h5ad",
  output_test_mod1 = "output/public_datasets/predict_modality/openproblems_bmmc_multiome_mod2/openproblems_bmmc_multiome_mod2.censor_dataset.output_test_mod1.h5ad",
  output_test_mod2 = "output/public_datasets/predict_modality/openproblems_bmmc_multiome_mod2/openproblems_bmmc_multiome_mod2.censor_dataset.output_test_mod2.h5ad",
  seed = 1L,
  max_mod1_columns = NULL,
  max_mod2_columns = 5000L
)
## VIASH END

cat("Reading input data\n")
ad1_raw <- anndata::read_h5ad(par$input_mod1)
ad2_raw <- anndata::read_h5ad(par$input_mod2)
ad1_mod <- unique(ad1_raw$var[["feature_types"]])
ad2_mod <- unique(ad2_raw$var[["feature_types"]])
new_dataset_id <- paste0(ad1_raw$uns[["dataset_id"]], "_PM_", tolower(ad1_mod), "2", tolower(ad2_mod))
common_uns <- list(dataset_id = new_dataset_id)

process_ad <- function(ad) {
  mod <- unique(ad$var[["feature_types"]])

  if (mod == "GEX") {
    cat("Applying precomputed size factors\n")
    ad_X <- as(ad$X, "CsparseMatrix")
    size_factors <-
      if (!is.null(ad$obs[["size_factors"]])) {
        ad$obs[["size_factors"]]
      } else {
        rep(1, nrow(ad))
      }
    ad_X@x <- log10(ad_X@x / size_factors[ad_X@i + 1] + 1)
    ad$layers <- list(counts = ad$X)
    ad$X <- ad_X
  } else if (mod == "ADT") {
    cat("CLR normalizing ADT data\n")
    obj <- Seurat::CreateSeuratObject(counts = ad$X, assay = "ADT")
    obj <- Seurat::NormalizeData(obj, normalization.method = "CLR", margin = 2)
    ad_X <- as(obj[["ADT"]]@data, "CsparseMatrix")
    ad$layers <- list(counts = ad$X)
    ad$X <- ad_X
  } else if (mod == "ATAC") {
    cat("Binarizing ATAC data\n")
    ad_X <- as(ad$X, "CsparseMatrix")
    ad_X@x <- (ad_X@x > 0) + 0
    ad$X <- ad_X
  }

  ad
}

ad1_raw <- process_ad(ad1_raw)
ad2_raw <- process_ad(ad2_raw)

# copied from scUtils/variance.R
colVars_spm <- function( spm ) {
  stopifnot( methods::is( spm, "dgCMatrix" ) )
  ans <- sapply( base::seq.int(spm@Dim[2]), function(j) {
    if( spm@p[j+1] == spm@p[j] ) { return(0) } # all entries are 0: var is 0
    mean <- base::sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) + mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) 
  })
  ans / ( spm@Dim[1] - 1 )
}

if (!is.null(par$max_mod1_columns) && par$max_mod1_columns < ncol(ad1_raw)) {
  cat("Sampling mod1 columns\n")
  ad1_var <- colVars_spm(ad1_raw$X)
  ad1_ix <- sample.int(ncol(ad1_raw), par$max_mod1_columns, prob = ad1_var)
  ad1_raw <- ad1_raw[, ad1_ix]
}
if (!is.null(par$max_mod2_columns) && par$max_mod2_columns < ncol(ad2_raw)) {
  cat("Sampling mod2 columns\n")
  ad2_var <- colVars_spm(ad2_raw$X)
  ad2_ix <- sample.int(ncol(ad2_raw), par$max_mod2_columns, prob = ad2_var)
  ad2_raw <- ad2_raw[, ad2_ix]
}

cat("Creating train objects\n")
ad1_var <- ad1_raw$var %>% select(one_of("gene_ids"), feature_types)
ad2_var <- ad2_raw$var %>% select(one_of("gene_ids"), feature_types)
is_train <- ad1_raw$obs[["is_train"]]
train_obs <- ad1_raw$obs[is_train, , drop = FALSE] %>% select(one_of("batch"))
test_obs <- ad1_raw$obs[!is_train, , drop = FALSE]  %>% select(one_of("batch"))

out_train_mod1 <- anndata::AnnData(
  X = ad1_raw$X[is_train, , drop = FALSE],
  layers = list(counts = ad1_raw$layers$counts[is_train, , drop = FALSE]),
  obs = train_obs,
  var = ad1_var,
  uns = common_uns
)
out_train_mod2 <- anndata::AnnData(
  X = ad2_raw$X[is_train, , drop = FALSE],
  layers = list(counts = ad2_raw$layers$counts[is_train, , drop = FALSE]),
  obs = train_obs,
  var = ad2_var,
  uns = common_uns
)

cat("Create test objects\n")
out_test_mod1 <- anndata::AnnData(
  X = ad1_raw$X[!is_train, , drop = FALSE],
  layers = list(counts = ad1_raw$layers$counts[!is_train, , drop = FALSE]),
  obs = test_obs,
  var = ad1_var,
  uns = common_uns
)
out_test_mod2 <- anndata::AnnData(
  X = ad2_raw$X[!is_train, , drop = FALSE],
  layers = list(counts = ad2_raw$layers$counts[!is_train, , drop = FALSE]),
  obs = test_obs,
  var = ad2_var,
  uns = common_uns
)

cat("Saving output files as h5ad\n")
zzz <- out_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- out_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- out_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- out_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")

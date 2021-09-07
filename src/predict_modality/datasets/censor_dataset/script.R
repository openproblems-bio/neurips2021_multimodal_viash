cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
# data_path <- "output/public_datasets/common/openproblems_bmmc_multiome/openproblems_bmmc_multiome.manual_formatting."
data_path <- "output/public_datasets/common/dyngen_citeseq_1/dyngen_citeseq_1.split_traintest."
out_path <- data_path %>% gsub("/common/", "/predict_modality/", .) %>% gsub("[^\\.]*\\.$", "censor_dataset\\.", .)
par <- list(
  input_mod1 = paste0(data_path, "output_mod2.h5ad"),
  input_mod2 = paste0(data_path, "output_rna.h5ad"),
  output_train_mod1 = paste0(out_path, "output_train_mod1.h5ad"),
  output_train_mod2 = paste0(out_path, "output_train_mod2.h5ad"),
  output_test_mod1 = paste0(out_path, "output_test_mod1.h5ad"),
  output_test_mod2 = paste0(out_path, "output_test_mod2.h5ad"),
  seed = 1L,
  max_mod1_columns = NULL,
  max_mod2_columns = 5000L
)
## VIASH END

cat("Reading input data\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)
ad1_mod <- unique(input_mod1$var[["feature_types"]])
ad2_mod <- unique(input_mod2$var[["feature_types"]])
new_dataset_id <- paste0(input_mod1$uns[["dataset_id"]], "_PM_", tolower(ad1_mod), "2", tolower(ad2_mod))
common_uns <- list(dataset_id = new_dataset_id)

process_ad <- function(ad) {
  mod <- unique(ad$var[["feature_types"]])
  ad$layers <- list(counts = as(ad$X, "CsparseMatrix"))
  
  # normalize
  if (mod == "GEX") {
    cat("Applying precomputed size factors\n")
    size_factors <-
      if (!is.null(ad$obs[["size_factors"]])) {
        ad$obs[["size_factors"]]
      } else {
        rep(1, nrow(ad))
      }
    ad$X@x <- log10(ad$X@x / size_factors[ad$X@i + 1] + 1)
  } else if (mod == "ADT") {
    cat("CLR normalizing ADT data\n")
    obj <- Seurat::CreateSeuratObject(counts = ad$X, assay = "ADT")
    obj <- Seurat::NormalizeData(obj, normalization.method = "CLR", margin = 2)
    ad$X <- as(obj[["ADT"]]@data, "CsparseMatrix")
  } else if (mod == "ATAC") {
    cat("Sampling and binarizing ATAC data\n")
    poss_ix <- which(colSums(ad$X) > 0)
    ad <- ad[, sort(sample(poss_ix, 10000))]$copy()
    ad$X@x <- (ad$X@x > 0) + 0
  }
  ad
}

input_mod1 <- process_ad(input_mod1)
input_mod2 <- process_ad(input_mod2)

cat("Creating train objects\n")
ad1_var <- input_mod1$var %>% select(one_of("gene_ids"), feature_types)
ad2_var <- input_mod2$var %>% select(one_of("gene_ids"), feature_types)
is_train <- input_mod1$obs[["is_train"]]
train_obs <- input_mod1$obs[is_train, , drop = FALSE] %>% select(one_of("batch"))
test_obs <- input_mod1$obs[!is_train, , drop = FALSE]  %>% select(one_of("batch"))

output_train_mod1 <- anndata::AnnData(
  X = input_mod1$X[is_train, , drop = FALSE],
  layers = list(counts = input_mod1$layers[["counts"]][is_train, , drop = FALSE]),
  obs = train_obs,
  var = ad1_var,
  uns = common_uns
)
output_train_mod2 <- anndata::AnnData(
  X = input_mod2$X[is_train, , drop = FALSE],
  layers = list(counts = input_mod2$layers[["counts"]][is_train, , drop = FALSE]),
  obs = train_obs,
  var = ad2_var,
  uns = common_uns
)

cat("Create test objects\n")
output_test_mod1 <- anndata::AnnData(
  X = input_mod1$X[!is_train, , drop = FALSE],
  layers = list(counts = input_mod1$layers[["counts"]][!is_train, , drop = FALSE]),
  obs = test_obs,
  var = ad1_var,
  uns = common_uns
)
output_test_mod2 <- anndata::AnnData(
  X = input_mod2$X[!is_train, , drop = FALSE],
  layers = list(counts = input_mod2$layers[["counts"]][!is_train, , drop = FALSE]),
  obs = test_obs,
  var = ad2_var,
  uns = common_uns
)

cat("Saving output files as h5ad\n")
zzz <- output_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- output_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- output_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- output_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")

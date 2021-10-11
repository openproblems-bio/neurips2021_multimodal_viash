cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
data_path <- "output/datasets/common/openproblems_bmmc_multiome_phase1/openproblems_bmmc_multiome_phase1.manual_formatting."
# data_path <- "output/public_datasets/common/dyngen_citeseq_1/dyngen_citeseq_1.split_traintest."
out_path <- data_path %>% gsub("/common/", "/predict_modality/", .) %>% gsub("[^\\.]*\\.$", "censor_dataset\\.", .)
par <- list(
  input_mod1 = paste0(data_path, "output_rna.h5ad"),
  input_mod2 = paste0(data_path, "output_mod2.h5ad"),
  output_train_mod1 = paste0(out_path, "output_train_mod1.h5ad"),
  output_train_mod2 = paste0(out_path, "output_train_mod2.h5ad"),
  output_test_mod1 = paste0(out_path, "output_test_mod1.h5ad"),
  output_test_mod2 = paste0(out_path, "output_test_mod2.h5ad"),
  seed = 1L
)
## VIASH END

set.seed(par$seed)

cat("Reading input data\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)
ad1_mod <- unique(input_mod1$var[["feature_types"]])
ad2_mod <- unique(input_mod2$var[["feature_types"]])
new_dataset_id <- paste0(input_mod1$uns[["dataset_id"]], "_PM_", tolower(ad1_mod), "2", tolower(ad2_mod))
common_uns <- list(dataset_id = new_dataset_id)


if (ad2_mod == "ATAC") {
  if (ncol(input_mod2) > 10000) {
    poss_ix <- which(Matrix::colSums(input_mod2$X) > 0)
    input_mod2 <- input_mod2[, sort(sample(poss_ix, 10000))]$copy()
  }
  input_mod2$X@x <- (input_mod2$X@x > 0) + 0
}

cat("Creating train objects\n")
ad1_var <- input_mod1$var %>% select(one_of("gene_ids"), feature_types)
ad2_var <- input_mod2$var %>% select(one_of("gene_ids"), feature_types)
is_train <- which(input_mod1$obs[["is_train"]])
is_test <- which(!input_mod1$obs[["is_train"]])

# sample cells
if (length(is_test) > 1000) {
  ct <- input_mod1$obs[["cell_type"]][is_test] %>% as.character()
  ct_tab <- table(ct)
  ct_freq <- setNames(as.vector(ct_tab) / sum(ct_tab), names(ct_tab))
  is_test <- sample(is_test, 1000, prob = sqrt(1 / ct_freq[ct]))
  # table(input_mod1$obs[["cell_type"]][is_test])
}

train_obs <- input_mod1$obs[is_train, , drop = FALSE] %>% select(one_of("batch", "size_factors"))
test_obs <- input_mod1$obs[is_test, , drop = FALSE] %>% select(one_of("batch", "size_factors"))

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
  X = input_mod1$X[is_test, , drop = FALSE],
  layers = list(counts = input_mod1$layers[["counts"]][is_test, , drop = FALSE]),
  obs = test_obs,
  var = ad1_var,
  uns = common_uns
)
output_test_mod2 <- anndata::AnnData(
  X = input_mod2$X[is_test, , drop = FALSE],
  layers = list(counts = input_mod2$layers[["counts"]][is_test, , drop = FALSE]),
  obs = test_obs,
  var = ad2_var,
  uns = common_uns
)

cat("Saving output files as h5ad\n")
zzz <- output_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- output_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- output_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- output_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")

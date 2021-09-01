cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
input_path <- "resources_test/common/test_resource."
# input_path <- "output/public_datasets/common/openproblems_bmmc_multiome/openproblems_bmmc_multiome.manual_formatting."
output_path <- ""
# output_path <- "output/public_datasets/match_modality/openproblems_bmmc_multiome/openproblems_bmmc_multiome.censor_dataset."

par <- list(
  input_mod1 = paste0(input_path, "output_rna.h5ad"),
  input_mod2 = paste0(input_path, "output_mod2.h5ad"),
  output_train_mod1 = paste0(output_path, "output_train_mod1.h5ad"),
  output_train_mod2 = paste0(output_path, "output_train_mod2.h5ad"),
  output_train_sol = paste0(output_path, "output_train_sol.h5ad"),
  output_test_mod1 = paste0(output_path, "output_test_mod1.h5ad"),
  output_test_mod2 = paste0(output_path, "output_test_mod2.h5ad"),
  output_test_sol = paste0(output_path, "output_test_sol.h5ad"),
  seed = 1L,
  knn = 10L
)
## VIASH END

set.seed(par$seed)

cat("Reading input data\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)

new_dataset_id <- paste0(input_mod1$uns[["dataset_id"]], "_MM")
common_uns <- list(dataset_id = new_dataset_id)

cat("Shuffle train cells\n")
train_ix <- which(input_mod1$obs$is_train) %>% sort
train_mod2_ix <- sample.int(length(train_ix))

cat("Shuffle test cells\n")
test_ix <- which(!input_mod1$obs$is_train) %>% sort
test_mod2_ix <- sample.int(length(test_ix))


cat("Creating train objects\n")
mod1_var <- input_mod1$var %>% select(one_of("gene_ids", "feature_types"))
mod2_var <- input_mod2$var %>% select(one_of("gene_ids", "feature_types"))
train_obs <- input_mod1$obs[train_ix, , drop = FALSE] %>% select(one_of("batch", "size_factors"))

output_train_mod1 <- anndata::AnnData(
  X = input_mod1$X[train_ix, , drop = FALSE],
  obs = train_obs,
  var = mod1_var,
  uns = common_uns
)
output_train_mod2 <- anndata::AnnData(
  X = input_mod2$X[train_ix[train_mod2_ix], , drop = FALSE] %>%
    magrittr::set_rownames(., paste0("cell_", seq_len(nrow(.)))),
  var = mod2_var,
  uns = common_uns
)

cat("Create test objects\n")
test_obs <- input_mod1$obs[test_ix, , drop = FALSE] %>% select(one_of("batch", "size_factors"))

output_test_mod1 <- anndata::AnnData(
  X = input_mod1$X[test_ix, , drop = FALSE],
  obs = test_obs,
  var = mod1_var,
  uns = common_uns
)
output_test_mod2 <- anndata::AnnData(
  X = input_mod2$X[test_ix[test_mod2_ix], , drop = FALSE] %>% 
    magrittr::set_rownames(., paste0("cell_", seq_len(nrow(.)))),
  var = mod2_var,
  uns = common_uns
)

cat("Create solution objects\n")
comp_neighbor_mat <- function(X1, X2) {
  dr1 <- lmds::lmds(X1, ndim = 100, distance_method = "spearman")
  dr2 <- lmds::lmds(X2, ndim = 100, distance_method = "spearman")
  knnidx1 <- cbind(seq_len(nrow(dr1)), FNN::knn.index(dr1, k = par$knn - 1))
  knncor1 <- do.call(rbind, map(seq_len(nrow(knnidx1)), function(i) {
    cor <- dynutils::calculate_similarity(
      x = X1[i, , drop = FALSE],
      y = X1[knnidx1[i, -1, drop = TRUE], , drop = FALSE],
      method = "spearman"
    )[1,]
    c(1, cor)
  }))
  knnidx2 <- cbind(seq_len(nrow(dr2)), FNN::knn.index(dr2, k = par$knn - 1))
  knncor2 <- do.call(rbind, map(seq_len(nrow(knnidx2)), function(i) {
    cor <- dynutils::calculate_similarity(
      x = X2[i, , drop = FALSE],
      y = X2[knnidx2[i, -1, drop = TRUE], , drop = FALSE],
      method = "spearman"
    )[1,]
    c(1, cor)
  }))

  bind_rows(
    tibble(i = as.vector(row(knnidx1)), j = as.vector(knnidx1), x = as.vector(knncor1)),
    tibble(i = as.vector(knnidx1),j = as.vector(row(knnidx1)), x = as.vector(knncor1)),
    tibble(i = as.vector(row(knnidx2)), j = as.vector(knnidx2), x = as.vector(knncor2)),
    tibble(i = as.vector(knnidx2), j = as.vector(row(knnidx2)), x = as.vector(knncor2))
  ) %>%
  group_by(i, j) %>%
  summarise(x = mean(x), .groups = "drop") %>% {
    Matrix::sparseMatrix(
      i = .$i, j = .$j, x = .$x, dims = list(nrow(dr1), nrow(dr1))
    )
  }
}

train_sol_mat <- Matrix::sparseMatrix(
  i = seq_along(train_mod2_ix),
  j = order(train_mod2_ix),
  x = rep(1, length(train_mod2_ix))
)
train_solknn <- comp_neighbor_mat(output_train_mod1$X, output_train_mod2$X)[, train_mod2_ix]
output_train_sol <- anndata::AnnData(
  X = train_sol_mat,
  obs = input_mod1$obs[train_ix, , drop = FALSE],
  layers = list(neighbors = train_solknn),
  uns = list(dataset_id = new_dataset_id, pairing_ix = train_mod2_ix - 1)
)

test_sol_mat <- Matrix::sparseMatrix(
  i = seq_along(test_mod2_ix),
  j = order(test_mod2_ix),
  x = rep(1, length(test_mod2_ix))
)
test_solknn <- comp_neighbor_mat(output_test_mod1$X, output_test_mod2$X)[, test_mod2_ix]
output_test_sol <- anndata::AnnData(
  X = test_sol_mat,
  obs = input_mod1$obs[test_ix, , drop = FALSE],
  layers = list(neighbors = test_solknn),
  uns = list(dataset_id = new_dataset_id, pairing_ix = test_mod2_ix - 1)
)

# checks
mean(rowSums(train_solknn > 0))
mean(rowSums(test_solknn > 0))
sum(train_solknn * train_sol_mat) == nrow(train_sol_mat)
sum(test_solknn * test_sol_mat) == nrow(test_sol_mat)

cat("Saving output files as h5ad\n")
zzz <- output_train_mod1$write_h5ad(par$output_train_mod1, compression = "gzip")
zzz <- output_train_mod2$write_h5ad(par$output_train_mod2, compression = "gzip")
zzz <- output_train_sol$write_h5ad(par$output_train_sol, compression = "gzip")
zzz <- output_test_mod1$write_h5ad(par$output_test_mod1, compression = "gzip")
zzz <- output_test_mod2$write_h5ad(par$output_test_mod2, compression = "gzip")
zzz <- output_test_sol$write_h5ad(par$output_test_sol, compression = "gzip")

cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = "resources_test/task1/pbmc_1k_protein_v3.solution.h5ad",
  input_prediction = "resources_test/task1/pbmc_1k_protein_v3.prediction.h5ad",
  output = "resources_test/task1/pbmc_1k_protein_v3.scores.h5ad"
)
## VIASH END

cat("Read solution h5ad\n")
adata_solution <- anndata::read_h5ad(par$input_solution)

cat("Read prediction h5ad\n")
expect_true(
  grepl("\\.h5ad$", par$input_prediction),
  info = "Prediction file should be an h5ad file"
)
adata_prediction <-
  tryCatch({
    anndata::read_h5ad(par$input_prediction)
  }, error = function(e) {
    stop(paste0("Can't open prediction h5ad file. Detailed error message:\n", e$message))
  })
expect_true(
  all.equal(dim(adata_solution), dim(adata_prediction)),
  info = "Dataset and prediction anndata objects should have the same shape / dimensions."
)
expect_true(
  adata_solution$uns$dataset_id == adata_prediction$uns$dataset_id
)

cat("Data wrangling\n")
tv <- as.matrix(adata_solution$X) %>% as.vector()
pv <- as.matrix(adata_prediction$X) %>% as.vector()

# visual check:
# qplot(tv, pv)

cat("Computing metrics\n")
rmse <- sqrt(mean((tv - pv) ** 2))
cor_pearson <- cor(tv, pv, method = "pearson")
cor_spearman <- cor(tv, pv, method = "spearman")

cat("Create output object\n")
out <- anndata::AnnData(
  X = NULL,
  shape = dim(adata_solution),
  uns = list(
    dataset_id = adata_prediction$uns$dataset_id,
    method_id = adata_prediction$uns$method_id,
    metric_ids = c("rmse", "cor_pearson", "cor_spearman"),
    metric_values = c(rmse, cor_pearson, cor_spearman),
    metric_moreisbetter = c(FALSE, TRUE, TRUE)
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

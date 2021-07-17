cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
# par <- list(
#   input_solution = c("resources_test/task1/test_resource.solution.h5ad"),
#   input_prediction = c("resources_test/task1/test_resource.prediction.h5ad"),
#   output = "test_resource.scores.h5ad"
# )
par <- list(
  input_solution = list.files("/tmp/neurips2021_work/d7/56ea677b97f0d7b03286e028749113", full.names = TRUE, recursive = TRUE, pattern = "*.output_solution.h5ad"),
  input_prediction = list.files("/tmp/neurips2021_work/d7/56ea677b97f0d7b03286e028749113", full.names = TRUE, recursive = TRUE, pattern = "*.output.h5ad"),
  output = "test_resource.scores.h5ad"
)
## VIASH END

cat("Reading solution info\n")
sol_info <- map_df(par$input_solution, function(path) {
  anndata::read_h5ad(path, backed = TRUE)$uns[c("dataset_id")] %>%
    as_tibble() %>% 
    mutate(sol_path = path)
})
cat("Reading prediction info\n")
pred_info <- map_df(par$input_prediction, function(path) {
  anndata::read_h5ad(path, backed = TRUE)$uns[c("dataset_id", "method_id")] %>%
    as_tibble() %>% 
    mutate(pred_path = path)
})
meta_info <- left_join(
  sol_info %>% crossing(method_id = unique(pred_info$method_id)),
  pred_info,
  by = c("dataset_id", "method_id")
)


cat("Evaluating predictions\n")
# list2env(dynutils::extract_row_to_list(meta_info, 1), .GlobalEnv)
out <- pmap_df(meta_info, function(dataset_id, method_id, sol_path, pred_path) {
  
  if (!is.na(pred_path)) {
    cat("Reading ", basename(sol_path), " and ", basename(pred_path), "\n", sep = "")
    # Read solution h5ad
    adata_solution <- anndata::read_h5ad(sol_path)

    # Read prediction h5ad
    expect_true(
      grepl("\\.h5ad$", pred_path),
      info = "Prediction file should be an h5ad file"
    )
    adata_prediction <-
      tryCatch({
        anndata::read_h5ad(pred_path)
      }, error = function(e) {
        stop(paste0("Can't open prediction h5ad file. Detailed error message:\n", e$message))
      })
    expect_true(
      all.equal(dim(adata_solution), dim(adata_prediction)),
      info = "Dataset and prediction anndata objects should have the same shape / dimensions."
    )

    # Wrangle data
    tv <- adata_solution$X
    pv <- adata_prediction$X

    # Compute metrics
    rmse <- sqrt(mean((tv - pv) ^ 2))
    score_pearson <- .5 - mean(diag(dynutils::calculate_similarity(tv, pv, method = "pearson", margin = 2, diag = TRUE, drop0 = TRUE))) / 2
    score_spearman <- .5 - mean(diag(dynutils::calculate_similarity(tv, pv, method = "spearman", margin = 2, diag = TRUE, drop0 = TRUE))) / 2
  } else {
    rmse <- Inf
    score_pearson <- 0
    score_spearman <- 0
  }

  tibble(
    dataset_id,
    method_id,
    metric_ids = c("rmse", "score_pearson", "score_spearman"),
    metric_values = c(rmse, score_pearson, score_spearman),
    metric_moreisbetter = c(FALSE, TRUE, TRUE)
  )
})


cat("Create output object\n")
out <- anndata::AnnData(
  X = NULL,
  shape = list(1L, 1L),
  uns = list(
    dataset_id = out$dataset_id,
    method_id = out$method_id,
    metric_ids = out$metric_ids,
    metric_values = out$metric_values,
    metric_moreisbetter = out$metric_moreisbetter
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

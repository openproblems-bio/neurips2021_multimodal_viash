## VIASH START
par <- list(
  input_solution = "resources/test/dyngen_bifurcating_antibody/dataset_task1_solution.h5ad",
  input_prediction = "resources/test/dyngen_bifurcating_antibody/dataset_task1_prediction_randomforest.h5ad",
  output = "resources/test/dyngen_bifurcating_antibody/dataset_task1_prediction_randomforest_scores.h5ad"
)
## VIASH END

###############################################################################
###                            LOAD DEPENDENCIES                            ###
###############################################################################

# load libraries
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(assertthat, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)


###############################################################################
###                             READ INPUT DATA                             ###
###############################################################################

# load dataset h5ad file
adata_solution <- anndata::read_h5ad(par$input_solution)

# load prediction h5ad file
assert_that(
  grepl("\\.h5ad$", par$input_prediction),
  msg = "Prediction file should be an h5ad file"
)
adata_prediction <-
  tryCatch({
    anndata::read_h5ad(par$input_prediction)
  }, error = function(e) {
    stop(paste0("Can't open prediction h5ad file. Detailed error message:\n", e$message))
  })
assert_that(
  all.equal(dim(adata_solution), dim(adata_prediction)),
  msg = "Dataset and prediction anndata objects should have the same shape / dimensions."
)
assert_that(
  adata_solution$uns$dataset_id == adata_prediction$uns$dataset_id
)

###############################################################################
###                           FETCH RELEVANT DATA                           ###
###############################################################################

# check input modalities
true_values <- adata_solution$layers["modality2"][adata_solution$obs$experiment == "test", adata_solution$var$modality2_has_values]

# check prediction layers
assert_that(
  adata_prediction$layers %has_name% "prediction",
  msg = "Prediction file needs to contain layer with name 'prediction'."
)
pred_values <- adata_prediction$layers[["prediction"]][adata_solution$obs$experiment == "test", adata_solution$var$modality2_has_values]

# calculate metrics
tv <- as.matrix(true_values) %>% as.vector()
pv <- as.matrix(pred_values) %>% as.vector()

# visual check:
# qplot(tv, pv)

###############################################################################
###                             COMPUTE METRICS                             ###
###############################################################################
rmse <- sqrt(mean((tv - pv) ** 2))
cor_pearson <- cor(tv, pv, method = "pearson")
cor_spearman <- cor(tv, pv, method = "spearman")

# tibble(
#   dataset_id = adata_prediction$uns$dataset_id,
#   method_id = adata_prediction$uns$method_id,
#   metric_id = c("rmse", "cor_pearson", "cor_spearman"),
#   metric_value = c(rmse, cor_pearson, cor_spearman)
# )
# return output file
out <- anndata::AnnData(
  X = NULL,
  shape = dim(adata_solution),
  uns = list(
    dataset_id = adata_prediction$uns$dataset_id,
    method_id = adata_prediction$uns$method_id,
    metric_id = c("rmse", "cor_pearson", "cor_spearman"),
    metric_value = c(rmse, cor_pearson, cor_spearman)
  )
)

out$write_h5ad(par$output, compression = "gzip")

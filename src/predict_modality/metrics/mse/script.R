cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
par <- list(
  input_solution = c("resources_test/predict_modality/test_resource.test_mod2.h5ad"),
  input_prediction = c("resources_test/predict_modality/test_resource.prediction.h5ad"),
  output = "test_resource.scores.h5ad"
)
## VIASH END

cat("Reading solution file\n")
ad_sol <- anndata::read_h5ad(par$input_solution)

cat("Reading prediction file\n")
ad_pred <- anndata::read_h5ad(par$input_prediction)

cat("Check prediction format\n")
expect_equal(
  ad_sol$uns$dataset_id, ad_pred$uns$dataset_id,
  info = "Prediction and solution have differing dataset_ids"
)
expect_true(
  all.equal(dim(ad_sol), dim(ad_pred)),
  info = "Dataset and prediction anndata objects should have the same shape / dimensions."
)

colVars_spm <- function( spm ) {
  stopifnot( methods::is( spm, "dgCMatrix" ) )
  ans <- sapply( base::seq.int(spm@Dim[2]), function(j) {
    if( spm@p[j+1] == spm@p[j] ) { return(0) } # all entries are 0: var is 0
    mean <- base::sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) + mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) 
  })
  ans / ( spm@Dim[1] - 1 )
}

cat("Computing MSE and MSLE metrics\n")
# Wrangle data
scores <- 
  full_join(
    summary(ad_sol$X) %>% rename(solx = x),
    summary(ad_pred$X) %>% rename(predx = x),
    by = c("i", "j")
  ) %>%
  select(-i, -j) %>%
  mutate(
    solx = ifelse(is.na(solx), 0, solx),
    predx = ifelse(is.na(predx), 0, predx),
  ) %>%
  summarise(
    mse = mean((solx - predx)^2)
  )

cat("Create output object\n")
out <- anndata::AnnData(
  uns = list(
    dataset_id = ad_pred$uns$dataset_id,
    method_id = ad_pred$uns$method_id,
    metric_ids = names(scores),
    metric_values = as.vector(as.matrix(scores))
  )
)

cat("Write output to h5ad file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

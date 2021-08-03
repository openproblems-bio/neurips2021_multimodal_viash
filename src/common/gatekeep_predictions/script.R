cat("Load dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(Matrix, quietly = TRUE, warn.conflicts = FALSE)
requireNamespace("anndata", quietly = TRUE)

## VIASH START
# par <- list(
#   input_solution = c("resources_test/predict_modality/test_resource.solution.h5ad"),
#   input_prediction = c("resources_test/predict_modality/test_resource.prediction.h5ad"),
#   output = "pairing.tsv"
# )
par <- list(
  input_solution = list.files("resources_test/joint_embedding", full.names = TRUE, recursive = TRUE, pattern = "*.solution.h5ad"),
  input_prediction = list.files("resources_test/joint_embedding", full.names = TRUE, recursive = TRUE, pattern = "*.prediction.h5ad"),
  output = "pairing.tsv"
)
## VIASH END

cat("Reading solution info\n")
sol_info <- map_df(par$input_solution, function(path) {
  tryCatch({
    anndata::read_h5ad(path, backed = TRUE)$uns[c("dataset_id")] %>%
      as_tibble() %>% 
      mutate(sol_path = path)
  }, error = function(e) {
    cat("Warning! Couldn't read from '", path, "'. Reason: \n", sep = "")
    cat(e$message)
    NULL
  })
})
cat("Reading prediction info\n")
pred_info <- map_df(par$input_prediction, function(path) {
  tryCatch({
    anndata::read_h5ad(path, backed = TRUE)$uns[c("dataset_id", "method_id")] %>%
      as_tibble() %>% 
      mutate(pred_path = path)
  }, error = function(e) {
    cat("Warning! Couldn't read from '", path, "'. Reason: \n", sep = "")
    cat(e$message)
    NULL
  })
})

cat("Join datasets with predictions\n")
meta_info <- left_join(
  sol_info %>% crossing(method_id = unique(pred_info$method_id)),
  pred_info,
  by = c("dataset_id", "method_id")
)

cat("\nPairing:\n")
print(meta_info)

unmatched1 <- anti_join(
  pred_info,
  sol_info %>% crossing(method_id = unique(pred_info$method_id)),
  by = c("dataset_id", "method_id")
)
if (nrow(unmatched1) > 0) {
  cat("\nWarning: Predictions found but no solutions:\n")
  print(unmatched1)
}

unmatched2 <- anti_join(
  sol_info %>% crossing(method_id = unique(pred_info$method_id)),
  pred_info,
  by = c("dataset_id", "method_id")
)
if (nrow(unmatched2) > 0) {
  cat("\nWarning: Solutions found but no predictions:\n")
  print(unmatched2)
}

cat("Write pairing info to file\n")
write_tsv(meta_info, par$output)

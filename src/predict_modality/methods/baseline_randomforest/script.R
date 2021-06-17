## VIASH START
par <- list(
  input = "resources/test/dyngen_bifurcating_antibody/dataset_task1_censor.h5ad",
  output = "output.h5ad"
)
## VIASH END

# load libraries
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)

library(randomForest, warn.conflicts = FALSE, quietly = TRUE)

# load h5ad file
adata <- anndata::read_h5ad(par$input)

# gather data
var_data <- adata$var
modality1 <- adata$layers[["modality1"]]
modality2 <- adata$layers[["modality2"]]

# perform DR on the features
dr <- prcomp(t(modality1), rank. = 10)$x

predictors1 <- t(modality1[,var_data$is_predictor])
responses1 <- t(modality1[,var_data$is_response])
predictors2 <- t(modality2[,var_data$is_predictor])

# predict for each gene
preds <- lapply(seq_len(ncol(predictors1)), function(i) {
  x <- cbind(dr[var_data$is_predictor,], expr = predictors1[,i])
  y <- predictors2[,1]

  rf <- randomForest::randomForest(x = x, y = y)

  newx <- cbind(dr[var_data$is_response,], expr = responses1[,i])
  predy <- stats::predict(rf, newx)

  predy
})

# create output
# TODO: create from scratch?
prediction <- Matrix::drop0(modality2 * 0)
prediction[, var_data$is_response] <- do.call(rbind, preds)

adata$layers[["prediction"]] <- prediction
adata$uns[["method_id"]] <- "baseline_randomforest"

adata$write_h5ad(par$output, compression = "gzip")

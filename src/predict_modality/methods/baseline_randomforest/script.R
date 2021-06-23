## VIASH START
par <- list(
  input = "resources/test/dyngen_bifurcating_antibody/dataset_task1_censor.h5ad",
  output = "resources/test/dyngen_bifurcating_antibody/dataset_task1_censor_randomforest.h5ad"
)
## VIASH END

# load libraries
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)
requireNamespace("randomForest", quietly = TRUE)

# load h5ad file
adata <- anndata::read_h5ad(par$input)

# gather data
obs_data <- adata$obs
modality1 <- adata$layers[["modality1"]]
modality2 <- adata$layers[["modality2"]]

# perform DR on the features
dr <- prcomp(modality1, rank. = 4)$x

dr_train <- dr[obs_data$experiment == "train",]
predictors_train <- modality1[obs_data$experiment == "train",]
responses_train <- modality2[obs_data$experiment == "train",]

dr_test <- dr[obs_data$experiment == "test",]
predictors_test <- modality1[obs_data$experiment == "test",]

# predict for each gene
preds <- lapply(seq_len(ncol(predictors_train)), function(i) {
  x <- cbind(dr_train, expr = predictors_train[,i])
  y <- responses_train[,i]

  if (length(unique(y)) > 1) {
    rf <- randomForest::randomForest(x = x, y = y)

    newx <- cbind(dr_test, expr = predictors_test[,i])
    stats::predict(rf, newx)
  } else {
    setNames(rep(unique(y), nrow(dr_test)), rownames(dr_test))
  }
})

# create output

prediction <- modality2 * 0
prediction[obs_data$experiment == "test"] <- do.call(cbind, preds)
prediction <- Matrix::drop0(prediction)

out <- anndata::AnnData(
  X = NULL,
  shape = dim(prediction),
  layers = list(
    prediction = prediction
  ),
  obs = adata$obs,
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    method_id = "baseline_randomforest"
  )
)

out$write_h5ad(par$output, compression = "gzip")

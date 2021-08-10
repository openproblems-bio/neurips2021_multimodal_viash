cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(keras, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_train_mod1 = "resources_test/match_modality/test_resource.train_mod1.h5ad",
  input_train_mod2 = "resources_test/match_modality/test_resource.train_mod2.h5ad",
  input_train_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  input_test_mod1 = "resources_test/match_modality/test_resource.test_mod1.h5ad",
  input_test_mod2 = "resources_test/match_modality/test_resource.test_mod2.h5ad",
  output = "output.h5ad",
  n_dims = 10L,
  distance_method = "spearman"
)
## VIASH END

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

match_train <- apply(input_train_sol$X, 1, function(x) which(x > 0))


cat("Running LMDS on input data\n")
# merge input matrices
mod1_X <- rbind(input_train_mod1$X, input_test_mod1$X)
mod2_X <- rbind(input_train_mod2$X[match_train, , drop = FALSE], input_test_mod2$X)

# perform DR
dr_x1 <- lmds::lmds(mod1_X, ndim = par$n_dims, distance_method = par$distance_method)
dr_x2 <- lmds::lmds(mod2_X, ndim = par$n_dims, distance_method = par$distance_method)

# split input matrices
dr_x1_train <- dr_x1[seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x2_train <- dr_x2[seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x1_test <- dr_x1[-seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x2_test <- dr_x2[-seq_len(nrow(input_train_mod1)), , drop = FALSE]


cat("Training neural network to predict mod2 dr\n")
model <-
  keras_model_sequential() %>%
  layer_dense(100, "relu", input_shape = par$n_dims) %>%
  layer_dense(32, "relu") %>%
  layer_dense(par$n_dims, "linear")

model %>% compile(
  loss = "mse",
  optimizer = "adam"
)

model %>% fit(dr_x1_train, dr_x2_train, epochs = 100, verbose = FALSE)
preds <- predict(model, dr_x1_test)
colnames(preds) <- colnames(dr_x2_test)


cat("Performing KNN between test mod2 DR and predicted test mod2\n")
par_frac <- 1
knn_out <- FNN::get.knnx(
  preds,
  dr_x2_test,
  k = min(100, ceiling(par_frac * nrow(preds)))
)

cat("Creating output data structures\n")
df <- tibble(
  i = as.vector(row(knn_out$nn.index)),
  j = as.vector(knn_out$nn.index),
  x = max(knn_out$nn.dist) * 2 - as.vector(knn_out$nn.dist)
)
knn_mat <- Matrix::sparseMatrix(i = df$i, j = df$j, x = df$x)

# # Plot mapping to see if this makes sense
# plotdf <- with(df %>% group_by(i) %>% slice(1), data.frame(dr_x2_test[i,], pred = preds[j,]))
# ggplot(plotdf) +
#   geom_point(aes(comp_1, comp_2, colour = "DR")) +
#   geom_point(aes(pred.comp_1, pred.comp_2, colour = "pred")) +
#   geom_segment(aes(x = comp_1, xend = pred.comp_1, y = comp_2, yend = pred.comp_2))


cat("Creating output anndata\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = input_train_mod1$uns[["dataset_id"]],
    method_id = "baseline_dr_nn_knn"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

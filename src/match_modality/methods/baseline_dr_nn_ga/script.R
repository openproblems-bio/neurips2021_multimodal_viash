cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("GA", quietly = TRUE)
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
  distance_method = "spearman",
  n_ga_pop = 200L,
  n_ga_iter = 500L
)
## VIASH END

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

match_train <- apply(input_train_sol$X, 1, function(x) which(x > 0)) %>% unname

cat("Running LMDS on input data\n")
# merge input matrices
mod1_X <- rbind(input_train_mod1$X, input_test_mod1$X)
mod2_X <- rbind(input_train_mod2$X[match_train, , drop = FALSE], input_test_mod2$X)

# perform DR
dr_x1 <- lmds::lmds(mod1_X, ndim = 10, distance_method = "pearson")
dr_x2 <- lmds::lmds(mod2_X, ndim = 3, distance_method = "pearson")

# split input matrices
dr_x1_train <- dr_x1[seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x2_train <- dr_x2[seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x1_test <- dr_x1[-seq_len(nrow(input_train_mod1)), , drop = FALSE]
dr_x2_test <- dr_x2[-seq_len(nrow(input_train_mod1)), , drop = FALSE]

cat("Training neural network to predict mod2 dr\n")
model <-
  keras_model_sequential() %>%
  layer_dense(100, "relu", input_shape = ncol(dr_x1)) %>%
  layer_dense(32, "relu") %>%
  layer_dense(ncol(dr_x2), "linear")

model %>% compile(
  loss = "mse",
  optimizer = "adam"
)

model %>% fit(dr_x1_train, dr_x2_train, epochs = 200, verbose = FALSE)
preds <- predict(model, dr_x1_test)
colnames(preds) <- colnames(dr_x2_test)



cat("Optimising correlation between mod1 distances and mod2 distances\n")
# x <- sample.int(nrow(dr_x2_test))
fitness <- function(x) {
  - mean((preds - dr_x2_test[x,])^2)
}
ga_out <- GA::ga(
  type = "permutation",
  fitness = fitness,
  lower = 1,
  upper = nrow(dr_x2_test),
  popSize = par$n_ga_pop,
  maxiter = par$n_ga_iter,
  # parallel = TRUE,
  monitor = GA:::gaMonitor
)

ord <- ga_out@solution[1,]
final_dr <- lmds::lmds(
  cbind(input_test_mod1$X, input_test_mod2$X[ord, , drop = FALSE]),
  ndim = par$n_dims,
  distance_method = par$distance_method
)
knn_out <- FNN::get.knn(final_dr, k = min(99, length(ord)-1))
knn_index <- cbind(seq_along(ord), knn_out$nn.index)
knn_dist <- cbind(rep(0, length(ord)), knn_out$nn.dist)

cat("Creating output data structures\n")
df <- tibble(
  i = as.vector(row(knn_out$nn.index)),
  j = as.vector(knn_out$nn.index),
  x = max(knn_out$nn.dist) * 2 - as.vector(knn_out$nn.dist)
)
knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  x = df$x,
  dims = list(nrow(dr_x1_test), nrow(dr_x2_test))
)

cat("Creating output anndata\n")
out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = input_train_mod1$uns[["dataset_id"]],
    method_id = method_id
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

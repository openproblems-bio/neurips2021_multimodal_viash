cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("GA", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(keras, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("batchelor", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
# path <- "resources_test/match_modality/test_resource."
# path <- "output/public_datasets/match_modality/dyngen_citeseq_1/dyngen_citeseq_1.censor_dataset.output_"
path <- "debug/debug."
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  input_train_sol = paste0(path, "train_sol.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_test_mod2 = paste0(path, "test_mod2.h5ad"),
  output = "output.h5ad",
  n_dims = 10L,
  distance_method = "spearman",
  n_ga_pop = 200L,
  n_ga_iter = 500L
)
meta <- list(functionality_name = "foo")

# # read in solution data to check whether method is working
# input_test_sol <- anndata::read_h5ad(paste0(path, "test_sol.h5ad"))
# match_test <- input_test_sol$uns$pairing_ix + 1
## VIASH END

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

cat("Log transform expression\n")
input_train_mod1$X@x <- log(input_train_mod1$X@x + 1)
input_train_mod2$X@x <- log(input_train_mod2$X@x + 1)
input_test_mod1$X@x <- log(input_test_mod1$X@x + 1)
input_test_mod2$X@x <- log(input_test_mod2$X@x + 1)

# reorder train_mod2 based on known solution
match_train <- input_train_sol$uns$pairing_ix + 1

# fetch batch labels
batch1 <- c(as.character(input_train_mod1$obs$batch), as.character(input_test_mod1$obs$batch))
# don't know batch ordering in input_test_mod2
batch2 <- c(as.character(input_train_mod1$obs$batch), rep("unknownbatch", nrow(input_test_mod2)))

cat("Running fastMNN\n")
mnn_out_1 <- batchelor::fastMNN(cbind(t(input_train_mod1$X), t(input_test_mod1$X)), batch = batch1)
dr_x1 <- SingleCellExperiment::reducedDim(mnn_out_1, "corrected")
mnn_out_2 <- batchelor::fastMNN(cbind(t(input_train_mod2$X[order(match_train), ]), t(input_test_mod2$X)), batch = batch2)
dr_x2 <- SingleCellExperiment::reducedDim(mnn_out_2, "corrected")
colnames(dr_x1) <- colnames(dr_x2) <- paste0("comp_", seq_len(ncol(dr_x1)))

# # visual checks
# ct1 <- c(as.character(input_train_sol$obs$cell_type), as.character(input_test_sol$obs$cell_type))
# ct2 <- c(as.character(input_train_sol$obs$cell_type), as.character(input_test_sol$obs$cell_type[match_test]))
# qplot(dr_x1[,1], dr_x1[,2], colour = factor(batch1))
# qplot(dr_x2[,1], dr_x2[,2], colour = factor(batch2))
# qplot(dr_x1[,1], dr_x1[,2], colour = factor(ct1))
# qplot(dr_x2[,1], dr_x2[,2], colour = factor(ct2))

# split DR matrices
train_ix <- seq_len(nrow(input_train_mod1))
dr_x1_train <- dr_x1[train_ix, , drop = FALSE]
dr_x1_test <- dr_x1[-train_ix, , drop = FALSE]
dr_x2_train <- dr_x2[train_ix, , drop = FALSE]
dr_x2_test <- dr_x2[-train_ix, , drop = FALSE]

cat("Training keras neural network to predict mod2 dr\n")
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

# # visual checks
# ggplot() +
#   geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(dr_x2_train, type = "train", input_train_sol$obs), size = 3) +
#   geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(dr_x2_test, type = "test", input_test_sol$obs[match_test, ]), size = 3) +
#   geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(preds, type = "pred", input_test_sol$obs), size = 3) +
#   facet_wrap(~type) +
#   theme_bw()


cat("Optimising correlation between mod1 distances and mod2 distances\n")
# x <- sample.int(nrow(dr_x2_test))
fitness <- function(x) {
  - mean((preds - dr_x2_test[x, ])^2)
}
ga_out <- GA::ga(
  type = "permutation",
  fitness = fitness,
  lower = 1,
  upper = nrow(dr_x2_test),
  popSize = par$n_ga_pop,
  maxiter = par$n_ga_iter,
  parallel = FALSE,
  monitor = GA:::gaMonitor,
  keepBest = TRUE
)
# ord <- ga_out@solution[1,]


bestSol <- do.call(rbind, ga_out@bestSol)
df <- reshape2::melt(bestSol, varnames = c("iter", "i"), value.name = "j") %>%
  group_by(i, j) %>%
  summarise(value = n(), .groups = "drop") %>%
  arrange(desc(value)) %>%
  head(input_test_mod1$n_obs * 1000)

# # visual checks
# ord <- order(match_test)
# plot_df <- data.frame(x = preds, dr_x2_test[ord,], input_test_sol$obs)
# ggplot(plot_df) + 
#   geom_segment(aes(x = comp_1, xend = x.comp_1, y = comp_2, yend = x.comp_2), alpha = .5) +
#   geom_point(aes(comp_1, comp_2, colour = cell_type, shape = "real"), size = 3) +
#   geom_point(aes(x.comp_1, x.comp_2, colour = cell_type, shape = "pred"), size = 3) +
#   theme_bw()

# final_dr <- lmds::lmds(
#   cbind(input_test_mod1$X, input_test_mod2$X[ord, , drop = FALSE]),
#   ndim = par$n_dims,
#   distance_method = par$distance_method
# )
# # ggplot() + 
# #   geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(final_dr, type = "real", input_test_sol$obs), size = 3) +
# #   theme_bw()
# knn_out <- FNN::get.knn(final_dr, k = min(999, length(ord)-1))
# knn_index <- cbind(seq_along(ord), knn_out$nn.index)
# knn_dist <- cbind(rep(0, length(ord)), knn_out$nn.dist)

# cat("Creating output data structures\n")
# df <- tibble(
#   i = as.vector(row(knn_index)),
#   j = as.vector(knn_index),
#   x = as.vector(knn_dist)
# ) %>% mutate(
#   j = order(ord)[j],
#   # rescale to get weights from distances
#   y = max(x) * 2 - x
# )
knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  # x = df$y,
  x = df$value,
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

cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_mod1 = "resources_test/match_modality/test_resource.mod1.h5ad",
  input_mod2 = "resources_test/match_modality/test_resource.mod2.h5ad",
  output = "output.h5ad",
  n_dims = 10L,
  n_waypoints = 100L,
  distance_method = "spearman",
  n_ga_pop = 50L,
  n_ga_iter = 20L
)
## VIASH END

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

# might want to perform lmds first?
# x1 <- ad1$X
# x2 <- ad2$X
x1 <- lmds::lmds(
  ad1$X,
  ndim = par$n_dims,
  distance_method = par$distance_method
)
x2 <- lmds::lmds(
  ad2$X,
  ndim = par$n_dims,
  distance_method = par$distance_method
)


waypoint_ix <- sample.int(nrow(x1), min(nrow(x1), par$n_waypoints))

cat("Calculating mod1 distances\n")
x1_dist <- dynutils::calculate_distance(
  x1[waypoint_ix, , drop = FALSE], 
  x1,
  method = "euclidean"
)

cat("Optimising correlation between mod1 distances and mod2 distances\n")
fitness <- function(x) {
  x2ord <- x2[x, ]
  x2_dist <- dynutils::calculate_distance(
    x2ord[waypoint_ix, , drop = FALSE], 
    x2ord,
    method = "euclidean"
  )
  cor(as.vector(x1_dist), as.vector(x2_dist), method = "pearson")
}
ga_out <- GA::ga(
  type = "permutation", 
  fitness = fitness,
  lower = 1,
  upper = nrow(x1),
  popSize = par$n_ga_pop,
  maxiter = par$n_ga_iter,
  #parallel = TRUE,
  monitor = GA:::gaMonitor
)

cat("Creating final mapping matrix\n")
ord <- ga_out@solution[1,]
final_dr <- lmds::lmds(
  cbind(ad1$X, ad2$X[ord, , drop = FALSE]),
  ndim = par$n_dims,
  distance_method = par$distance_method
)
knn_out <- FNN::get.knn(final_dr, k = 99)
knn_index <- cbind(seq_along(ord), knn_out$nn.index)
knn_dist <- cbind(rep(0, length(ord)), knn_out$nn.dist)

df <- tibble(
  i = as.vector(row(knn_index)),
  j = ord[as.vector(knn_index)],
  x = max(knn_dist) * 2 - as.vector(knn_dist)
)
knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  x = df$x
)

out <- anndata::AnnData(
  X = as(knn_mat, "CsparseMatrix"),
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = "baseline_optimize_distances"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

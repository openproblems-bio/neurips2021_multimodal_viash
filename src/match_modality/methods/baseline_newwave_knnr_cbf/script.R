cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("pbapply", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("NewWave", quietly = TRUE)
requireNamespace("FNN", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
# path <- "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# path <- "output/public_datasets/match_modality/dyngen_citeseq_1/dyngen_citeseq_1.censor_dataset.output_"
# path <- "output/public_datasets/match_modality/dyngen_atac_1/dyngen_atac_1.censor_dataset.output_"
# path <- "output/public_datasets/match_modality/dyngen_citeseq_3/dyngen_citeseq_3.censor_dataset.output_"
path <- "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# path <- "debug/debug."
par <- list(
  input_train_mod1 = paste0(path, "train_mod1.h5ad"),
  input_train_mod2 = paste0(path, "train_mod2.h5ad"),
  input_train_sol = paste0(path, "train_sol.h5ad"),
  input_test_mod1 = paste0(path, "test_mod1.h5ad"),
  input_test_mod2 = paste0(path, "test_mod2.h5ad"),
  output = "output.h5ad",
  n_pop = 300L,
  newwave_maxiter = 10
)
meta <- list(functionality_name = "foo")

# # read in solution data to check whether method is working
# input_test_sol <- anndata::read_h5ad(paste0(path, "test_sol.h5ad"))
# match_test <- input_test_sol$uns$pairing_ix + 1
## VIASH END

n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

# reorder train_mod2 based on known solution
match_train <- input_train_sol$uns$pairing_ix + 1

# fetch batch labels
batch1 <- c(as.character(input_train_mod1$obs$batch), as.character(input_test_mod1$obs$batch))
# don't know batch ordering in input_test_mod2
batch2 <- c(as.character(input_train_mod1$obs$batch), rep("unknownbatch", nrow(input_test_mod2)))

cat("Running NewWave\n")
data1 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = cbind(t(input_train_mod1$X), t(input_test_mod1$X))),
  colData = data.frame(batch = factor(batch1))
)
data1 <- data1[Matrix::rowSums(SummarizedExperiment::assay(data1)) > 0, ]
# option 1: filter by HVG
# data1 <- data1[order(proxyC::rowSds(SummarizedExperiment::assay(data1)), decreasing = TRUE)[1:100], ]

# option 2: use zinbsurf instead of zinbwave
# res1 <- zinbwave::zinbsurf(data1, K=10, X=~batch)
res1 <- NewWave::newWave(
  data1,
  X = "~batch",
  verbose = TRUE,
  K = 10,
  maxiter_optimize = par$newwave_maxiter,
  n_gene_par = min(par$newwave_ngene, nrow(data1)),
  n_cell_par = min(par$newwave_ncell, ncol(data1)),
  commondispersion = FALSE
)
dr_x1 <- SingleCellExperiment::reducedDim(res1)

data2 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = cbind(t(input_train_mod2$X[order(match_train), , drop = FALSE]), t(input_test_mod2$X))),
  colData = data.frame(batch = factor(batch2))
)
data2 <- data2[Matrix::rowSums(SummarizedExperiment::assay(data2)) > 0, ]
# data2 <- data2[order(proxyC::rowSds(SummarizedExperiment::assay(data2)), decreasing = TRUE)[1:100], ]

# res2 <- zinbwave::zinbsurf(data2, K=10, X=~batch)
res2 <- NewWave::newWave(
  data2,
  X = "~batch",
  verbose = TRUE,
  K = 10,
  maxiter_optimize = par$newwave_maxiter,
  n_gene_par = min(par$newwave_ngene, nrow(data2)),
  n_cell_par = min(par$newwave_ncell, ncol(data2)),
  commondispersion = FALSE
)
dr_x2 <- SingleCellExperiment::reducedDim(res2)

colnames(dr_x1) <- paste0("comp_", seq_len(ncol(dr_x1)))
colnames(dr_x2) <- paste0("comp_", seq_len(ncol(dr_x2)))

# # visual checks
# ct1 <- c(as.character(input_train_sol$obs$cell_type), as.character(input_test_sol$obs$cell_type))
# ct2 <- c(as.character(input_train_sol$obs$cell_type), as.character(input_test_sol$obs$cell_type[match_test]))
# patchwork::wrap_plots(
#   qplot(dr_x1[,1], dr_x1[,2], colour = factor(batch1)),
#   qplot(dr_x2[,1], dr_x2[,2], colour = factor(batch2)),
#   qplot(dr_x1[,1], dr_x1[,2], colour = factor(ct1)),
#   qplot(dr_x2[,1], dr_x2[,2], colour = factor(ct2))
# )


# split DR matrices
train_ix <- seq_len(nrow(input_train_mod1))
dr_x1_train <- dr_x1[train_ix, , drop = FALSE]
dr_x1_test <- dr_x1[-train_ix, , drop = FALSE]
dr_x2_train <- dr_x2[train_ix, , drop = FALSE]
dr_x2_test <- dr_x2[-train_ix, , drop = FALSE]

cat("Predicting mod1 DR of test cells\n")
pred_mod1 <- apply(dr_x1_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_x2_train,
    test = dr_x2_test,
    y = yi,
    k = min(15, nrow(dr_x1_test))
  )$pred
})
cat("Predicting mod2 DR of test cells\n")
pred_mod2 <- apply(dr_x2_train, 2, function(yi) {
  FNN::knn.reg(
    train = dr_x1_train,
    test = dr_x1_test,
    y = yi,
    k = min(15, nrow(dr_x1_test))
  )$pred
})


# # visual checks
# patchwork::wrap_plots(
#   ggplot() +
#     geom_point(aes(comp_1, comp_2, colour = cell_type), data.frame(rbind(dr_x2, pred_mod2)), size = 3, colour = "gray") +
#     geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(dr_x2_train, type = "1train", input_train_sol$obs), size = 3) +
#     geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(dr_x2_test, type = "2test", input_test_sol$obs[match_test, ]), size = 3) +
#     geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(pred_mod2, type = "3pred", input_test_sol$obs), size = 3) +
#     facet_wrap(~type) +
#     theme_bw(),
#   ggplot() +
#     geom_point(aes(comp_1, comp_2, colour = cell_type), data.frame(rbind(dr_x1, pred_mod1)), size = 3, colour = "gray") +
#     geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(dr_x1_train, type = "1train", input_train_sol$obs), size = 3) +
#     geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(dr_x1_test, type = "2test", input_test_sol$obs[match_test, ]), size = 3) +
#     geom_point(aes(comp_1, comp_2, colour = cell_type, shape = type), data.frame(pred_mod1, type = "3pred", input_test_sol$obs), size = 3) +
#     facet_wrap(~type) +
#     theme_bw(),
#   ncol = 1
# )


cat("Minimising distances between mod1 and mod2 pairs\n")
gen_vec <- function(z) {
  int <- seq_len(nrow(pred_mod1))

  i <- j <- c()
  resti <- int
  restj <- int

  while (length(resti) > 0) {
    ixi <- sample.int(length(resti), 1)
    newi <- resti[[ixi]]
    d1 <- proxy::dist(pred_mod1[restj, , drop = FALSE], dr_x1_test[newi, , drop = FALSE], method = "euclidean")
    d2 <- proxy::dist(pred_mod2[restj, , drop = FALSE], dr_x2_test[newi, , drop = FALSE], method = "euclidean")
    d12 <- d1 + d2
    ixj <- which.min(d12[, 1])
    newj <- restj[[ixj]]
    resti <- resti[-ixi]
    restj <- restj[-ixj]
    i <- c(i, newi)
    j <- c(j, newj)

    #  tibble(i, j); tibble(resti, restj)
  }

  tibble::tibble(i, j)
}

outs <- pbapply::pblapply(seq_len(par$n_pop), cl = n_cores, gen_vec)
# outs <- lapply(seq_len(par$n_pop), gen_vec)
df <- bind_rows(outs) %>%
  group_by(i, j) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  mutate(gold = i == j)

knn_mat <- Matrix::sparseMatrix(
  i = df$i,
  j = df$j,
  x = df$n,
  dims = list(nrow(dr_x1_test), nrow(dr_x2_test))
)

# normalise to make rows sum to 1
rs <- Matrix::rowSums(knn_mat)
knn_mat@x <- knn_mat@x / rs[knn_mat@i + 1]

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

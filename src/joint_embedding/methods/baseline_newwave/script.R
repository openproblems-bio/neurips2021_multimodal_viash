cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("NewWave", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
path <- "resources_test/joint_embedding/test_resource."
path <- "output/public_datasets/joint_embedding/dyngen_citeseq_1/dyngen_citeseq_1.censor_dataset.output_"
par <- list(
  input_mod1 = paste0(path, "mod1.h5ad"),
  input_mod2 = paste0(path, "mod2.h5ad"),
  output = "output.h5ad"
)
meta <- list(functionality_name = "foo")
# input_sol <- anndata::read_h5ad(paste0(path, "solution.h5ad"))
## VIASH END

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)
batch <- as.character(input_mod1$obs$batch)

cat("Running NewWave on mod1\n")
data1 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = t(input_mod1$X)),
  colData = data.frame(batch = factor(batch))
)
res1 <- NewWave::newWave(
  data1,
  X = "~batch",
  verbose = TRUE,
  K = 10,
  maxiter_optimize = 100,
  n_gene_par = min(300, nrow(data1)),
  n_cell_par = min(300, ncol(data1)),
  commondispersion = FALSE
)
dr_x1 <- SingleCellExperiment::reducedDim(res1)

cat("Running NewWave on mod2\n")
data2 <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = t(input_mod2$X)),
  colData = data.frame(batch = factor(batch))
)
res2 <- NewWave::newWave(
  data2,
  X = "~batch",
  verbose = TRUE,
  K = 10,
  maxiter_optimize = 100,
  n_gene_par = min(300, nrow(data2)),
  n_cell_par = min(300, ncol(data2)),
  commondispersion = FALSE
)
dr_x2 <- SingleCellExperiment::reducedDim(res2)

dr <- do.call(cbind, lapply(seq_len(ncol(dr_x1)), function(i) {
  cbind(dr_x1[,i], dr_x2[,i])
}))

# qplot(dr[,1], dr[,2], colour = batch)

rownames(dr) <- rownames(input_mod1)
colnames(dr) <- paste0("comp_", seq_len(ncol(dr)))

out <- anndata::AnnData(
  X = dr,
  uns = list(
    dataset_id = input_mod1$uns[["dataset_id"]],
    method_id = meta$functionality_name
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

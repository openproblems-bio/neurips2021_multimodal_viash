cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("batchelor", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
path <- "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter."
# path <- "output/public_datasets/joint_embedding/dyngen_citeseq_1/dyngen_citeseq_1.censor_dataset.output_"
par <- list(
  input_mod1 = paste0(path, "mod1.h5ad"),
  input_mod2 = paste0(path, "mod2.h5ad"),
  output = "output.h5ad"
)
meta <- list(functionality_name = "foo")
## VIASH END

method_id <- meta$functionality_name

cat("Reading h5ad files\n")
input_mod1 <- anndata::read_h5ad(par$input_mod1)
input_mod2 <- anndata::read_h5ad(par$input_mod2)

cat("Log transform expression\n")
input_mod1$X@x <- log(input_mod1$X@x + 1)
input_mod2$X@x <- log(input_mod2$X@x + 1)

cat("Running fastMNN\n")
mnn_out <- batchelor::fastMNN(
  rbind(t(input_mod1$X), t(input_mod2$X)), 
  batch = input_mod1$obs$batch
)
dr <- SingleCellExperiment::reducedDim(mnn_out, "corrected")

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

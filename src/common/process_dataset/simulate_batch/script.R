cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_rna = "resources_test/common/test_resource.output_rna.h5ad",
  input_mod2 = "resources_test/common/test_resource.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad",
  num_batches = 3L,
  batch_sd = .4
)
## VIASH END

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

if (!("batch" %in% colnames(ad_rna$obs) && "batch" %in% colnames(ad_mod2$obs))) {
  cat("Simulating batch effects\n")

  batch_names <- paste0("batch", seq_len(par$num_batches))
  batch <- sample(batch_names, nrow(ad_rna), replace = TRUE)

  weights_rna <- map(batch_names, function(name) {
    wt <- pmax(rnorm(ncol(ad_rna), mean = 1, sd = par$batch_sd), 0)
    Matrix::Matrix(diag(wt), sparse = TRUE)
  })
  weights_mod2 <- map(batch_names, function(name) {
    wt <- pmax(rnorm(ncol(ad_mod2), mean = 1, sd = par$batch_sd), 0)
    Matrix::Matrix(diag(wt), sparse = TRUE)
  })
  names(weights_rna) <- names(weights_mod2) <- batch_names

  simbatch_rna <- do.call(rbind, map(batch_names, function(bn) {
    (ad_rna$X[batch == bn, ] %*% weights_rna[[bn]]) %>%
      round %>%
      Matrix::drop0()
  }))
  simbatch_mod2 <- do.call(rbind, map(batch_names, function(bn) {
    (ad_mod2$X[batch == bn, ] %*% weights_mod2[[bn]]) %>%
      round %>%
      Matrix::drop0()
  }))

  ad_rna$X <- simbatch_rna[rownames(ad_rna), ]
  ad_mod2$X <- simbatch_mod2[rownames(ad_rna), ]
  ad_rna$obs[["batch"]] <- batch
  ad_mod2$obs[["batch"]] <- batch
}

# dr <- SCORPIUS::reduce_dimensionality(ad_rna$X)
# qplot(dr[,1], dr[,2], colour = ad_rna$obs[["batch"]])

cat("Writing mod1 data\n")
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

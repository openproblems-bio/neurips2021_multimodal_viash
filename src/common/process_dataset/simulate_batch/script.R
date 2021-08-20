cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input_rna = "resources_test/common/test_resource.output_rna.h5ad",
  input_mod2 = "resources_test/common/test_resource.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad",
  num_batches = 3L,
  batch_sd = .4
)
par <- list(
  input_rna = "work/0a/fec79afc2a21b0c6d1ef490d53e6be/totalvi_10x_malt_10k.pseudotime_order.output_rna.h5ad",
  input_mod2 = "work/0a/fec79afc2a21b0c6d1ef490d53e6be/totalvi_10x_malt_10k.pseudotime_order.output_mod2.h5ad",
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
  batch_ix <- sample.int(par$num_batches, nrow(ad_rna), replace = TRUE)

  # map rna counts
  summ_rna <- summary(ad_rna$X)
  weights_rna <- do.call(rbind, map(batch_names, function(name) {
    pmax(rnorm(ncol(ad_rna), mean = 1, sd = par$batch_sd), 0)
  }))
  summ_rna$bi <- batch_ix[summ_rna$i]
  summ_rna$mx <- round(summ_rna$x * weights_rna[cbind(summ_rna$bi, summ_rna$j)])
  ad_rna$X@x <- summ_rna$mx

  # map mod2 counts
  summ_mod2 <- summary(ad_mod2$X)
  weights_mod2 <- do.call(rbind, map(batch_names, function(name) {
    pmax(rnorm(ncol(ad_mod2), mean = 1, sd = par$batch_sd), 0)
  }))
  summ_mod2$bi <- batch_ix[summ_mod2$i]
  summ_mod2$mx <- round(summ_mod2$x * weights_mod2[cbind(summ_mod2$bi, summ_mod2$j)])
  ad_mod2$X@x <- summ_mod2$mx

  ad_rna$obs[["batch"]] <- batch_names[batch_ix]
  ad_mod2$obs[["batch"]] <- batch_names[batch_ix]
}

# dr <- SCORPIUS::reduce_dimensionality(ad_rna$X)
# qplot(dr[,1], dr[,2], colour = ad_rna$obs[["batch"]])

cat("Writing mod1 data\n")
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

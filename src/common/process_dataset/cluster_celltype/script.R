cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("lmds", quietly = TRUE)
requireNamespace("mclust", quietly = TRUE)

## VIASH START
par <- list(
  input_rna = "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_rna.h5ad",
  input_mod2 = "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad"
)
## VIASH END

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

if (!"cell_type" %in% colnames(ad_rna$obs)) {
  cat("Computing dimensionality reduction\n")
  dr <- lmds::lmds(cbind(ad_rna$X, ad_mod2$X), distance_method = "spearman")

  cat("Clustering with Mclust\n")
  mclustBIC <- mclust::mclustBIC
  labels <- paste0("cluster", mclust::Mclust(dr, verbose = TRUE)$classification)

  ad_rna$obs[["cell_type"]] <- labels
  ad_mod2$obs[["cell_type"]] <- labels
}

cat("Writing mod1 data\n")
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

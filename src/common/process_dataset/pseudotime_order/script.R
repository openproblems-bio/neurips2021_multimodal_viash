cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("SCORPIUS", quietly = TRUE)

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

obs_col_mod1 <- "pseudotime_order_GEX" # paste0("pseudotime_order_", unique(ad_mod1$var[["feature_types"]]))
obs_col_mod2 <- paste0("pseudotime_order_", unique(ad_mod2$var[["feature_types"]]))

if (!obs_col_mod1 %in% colnames(ad_rna$obs)) {
  cat("Computing pseudotime ordering for mod1\n")
  dr_rna <- SCORPIUS::reduce_dimensionality(as(ad_rna$X, "CsparseMatrix"))
  pt_rna <- SCORPIUS::infer_trajectory(dr_rna)$time
  ad_rna$obs[[obs_col_mod1]] <- pt_rna
  ad_mod2$obs[[obs_col_mod1]] <- pt_rna
}
if (!obs_col_mod2 %in% colnames(ad_mod2$obs)) {
  cat("Computing pseudotime ordering for mod2\n")
  dr_mod2 <- SCORPIUS::reduce_dimensionality(as(ad_mod2$X, "CsparseMatrix"))
  pt_mod2 <- SCORPIUS::infer_trajectory(dr_mod2)$time
  ad_rna$obs[[obs_col_mod2]] <- pt_mod2
  ad_mod2$obs[[obs_col_mod2]] <- pt_mod2
}

cat("Writing mod1 data\n")
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

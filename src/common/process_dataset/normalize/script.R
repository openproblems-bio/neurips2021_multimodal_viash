##
## TODO: Note that this script contains only very basic preprocessing.
##   This should be replaced at a later stage with something that makes more sense.
##

cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("scran", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)
sc <- reticulate::import("scanpy")

## VIASH START
par <- list(
  input_rna = "resources_test/common/test_resource.tmp2.output_rna.h5ad",
  input_mod2 = "resources_test/common/test_resource.tmp2.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad"
)
## VIASH END

cat("Reading mod1 h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)


cat("Computing clusters for scran\n")
ad_rna_pp <- ad_rna$copy()
sc$pp$normalize_per_cell(ad_rna_pp, counts_per_cell_after=1e6)
sc$pp$log1p(ad_rna_pp)
sc$pp$pca(ad_rna_pp, n_comps = 15L)
sc$pp$neighbors(ad_rna_pp)
sc$tl$leiden(ad_rna_pp, key_added = "groups", resolution = 0.5)
input_groups <- ad_rna_pp$obs$groups
rm(ad_rna_pp)

cat("Computing size factors\n")
mat <- t(ad_rna$X)
colnames(mat) <- as.vector(colnames(mat))
rownames(mat) <- as.vector(rownames(mat))
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
size_factors <- SingleCellExperiment::sizeFactors(scran::computeSumFactors(
  sce,
  clusters = input_groups, 
  min.mean = 0.1
))

ad_rna$obs[["size_factors"]] <- size_factors

cat("Reading mod2 h5ad file\n")
ad_mod2 <- anndata::read_h5ad(par$input_mod2)
mod2_type <- unique(ad_mod2$var$feature_types)

if (mod2_type == "ATAC") {
  cat("Binarizing peaks\n")
  ad_mod2$X@x <- (ad_mod2$X@x > 0) + 0
}


cat("Writing mod1 data\n")
print(ad_rna)
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
print(ad_mod2)
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

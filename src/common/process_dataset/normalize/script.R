##
## TODO: Note that this script contains only very basic preprocessing.
##   This should be replaced at a later stage with something that makes more sense.
##

cat("Loading dependencies\n")
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("scran", quietly = TRUE)
requireNamespace("SingleCellExperiment", quietly = TRUE)

## VIASH START
# par <- list(
#   input_rna = "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.tmp2.output_rna.h5ad",
#   input_mod2 = "resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.tmp2.output_mod2.h5ad",
#   output_rna = "output_rna.h5ad",
#   output_mod2 = "output_mod2.h5ad"
# )
# par <- list(
#   input_rna = "/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/work/c8/17106bcf362cbecf431db6f30bc475/totalvispleenlymph_spleen_lymph_206.quality_control.output_rna.h5ad",
#   input_mod2 = "/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/work/c8/17106bcf362cbecf431db6f30bc475/totalvispleenlymph_spleen_lymph_206.quality_control.output_mod2.h5ad",
#   output_rna = "output_rna.h5ad",
#   output_mod2 = "output_mod2.h5ad"
# )
par <- list(
  input_rna = "/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/work/e9/e66c4cdd574423ff9de2f3f825029b/10x_pbmc_1k_protein_v3.quality_control.output_rna.h5ad",
  input_mod2 = "/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/work/e9/e66c4cdd574423ff9de2f3f825029b/10x_pbmc_1k_protein_v3.quality_control.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad"
)
## VIASH END

cat("Reading mod1 h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)

cat("Log normalizing\n")
X <- ad_rna$X
X@x <- log(X@x + 1)

cat("Dimensionality reduction with LMDS\n")
dr <- lmds::lmds(X, distance_method = "euclidean", ndim = 15)

cat("Computing clusters for scran\n")
input_groups <- kmeans(dr, 15)$cluster

# cat("Computing clusters for scran\n")
# ad_rna_pp <- ad_rna$copy()
# cat("  Normalizing\n")
# zz <- sc$pp$normalize_per_cell(ad_rna_pp, counts_per_cell_after=1e6)
# cat("  Log1p\n")
# zz <- sc$pp$log1p(ad_rna_pp)
# cat("  PCA\n")
# # ad_rna_pp$obsm[["X_pca"]] <- lmds::lmds(ad_rna_pp$X, ndim = 15, distance_method = "pearson")
# zz <- sc$pp$pca(ad_rna_pp, n_comps = 15L)
# cat("  KNN\n")
# zz <- sc$pp$neighbors(ad_rna_pp)
# cat("  Leiden\n")
# input_groups <- leiden::leiden(as(ad_rna_pp$obsp[["connectivities"]], "CsparseMatrix"))
# # zz <- sc$tl$leiden(ad_rna_pp, key_added = "groups")
# # input_groups <- ad_rna_pp$obs$groups
# rm(ad_rna_pp)

cat("Computing size factors\n")
mat <- t(ad_rna$X)
colnames(mat) <- as.vector(colnames(mat))
rownames(mat) <- as.vector(rownames(mat))
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
size_factors <- SingleCellExperiment::sizeFactors(scran::computeSumFactors(
  sce,
  clusters = input_groups, 
  min.mean = 0.5
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

##
## TODO: Note that this script contains only very basic preprocessing.
##   This should be replaced at a later stage with something that makes more sense.
##

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
sc <- reticulate::import("scanpy")

## VIASH START
par <- list(
  input_rna = "output/pbmc_1k_protein_v3/pbmc_1k_protein_v3.download_10x_dataset.h5ad",
  input_mod2 = "output/pbmc_1k_protein_v3/pbmc_1k_protein_v3.download_10x_dataset.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad",
)
## VIASH END

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

# cat("Filtering cells\n")
# sc$pp$filter_cells(ad, min_counts = 100)
# 
cat("Filtering genes\n")
ad_rna$pp$filter_genes(ad_rna, min_counts = 100)
ad_mod2$pp$filter_genes(ad_mod2, min_counts = 100)
# 
# if ("chromatin" %in% ad$obsm_keys()) {
#   cat("Filtering chromatin values\n")
#   dat <- ad$obsm[["chromatin"]]
#   ix <- colSums(dat != 0) >= 10 # at least 10 non-zero cells
#   ad$obsm[["chromatin"]] <- as(dat[, ix, drop = FALSE], "RsparseMatrix")
#   ad$uns[["chromatin_varnames"]] <- ad$uns[["chromatin_varnames"]][ix]
# }
# if ("protein" %in% ad$obsm_keys()) {
#   cat("Filtering protein values\n")
#   dat <- ad$obsm[["protein"]]
#   ix <- colSums(dat != 0) >= 10 # at least 10 non-zero cells
#   ad$obsm[["protein"]] <- as(dat[, ix, drop = FALSE], "RsparseMatrix")
#   ad$uns[["protein_varnames"]] <- ad$uns[["protein_varnames"]][ix]
# }

cat("Storing output to '", par$output_rna, "'\n", sep = "")
ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Storing output to '", par$output_mod2, "'\n", sep = "")
ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

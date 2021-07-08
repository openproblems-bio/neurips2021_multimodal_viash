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
  input_rna = "output/common_datasets/pbmc_1k_protein_v3/pbmc_1k_protein_v3.output_rna.h5ad",
  input_mod2 = "output/common_datasets/pbmc_1k_protein_v3/pbmc_1k_protein_v3.output_mod2.h5ad",
  output_rna = "output/common_datasets/pbmc_1k_protein_v3/pbmc_1k_protein_v3.normalized.output_rna.h5ad",
  output_mod2 = "output/common_datasets/pbmc_1k_protein_v3/pbmc_1k_protein_v3.normalized.output_mod2.h5ad"
)
## VIASH END

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

cat("Filtering genes\n")
sc$pp$filter_genes(ad_rna, min_counts = 100)
sc$pp$filter_genes(ad_mod2, min_counts = 100)

cat("Filtering cells\n")
# sc$pp$filter_cells(ad, min_counts = 100)
mat <- cbind(ad_rna$X, ad_mod2$X)
filt <- rowSums(mat) > 100
ad_rna <- ad_rna[filt, ]
ad_mod2 <- ad_mod2[filt, ]

cat("Storing output to '", par$output_rna, "'\n", sep = "")
ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Storing output to '", par$output_mod2, "'\n", sep = "")
ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

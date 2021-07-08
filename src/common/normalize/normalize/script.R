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
  input_rna = "resources_test/common/pbmc_1k_protein_v3.output_rna.h5ad",
  input_mod2 = "resources_test/common/pbmc_1k_protein_v3.output_mod2.h5ad",
  output_rna = "resources_test/commons/pbmc_1k_protein_v3.normalize.output_rna.h5ad",
  output_mod2 = "resources_test/common/pbmc_1k_protein_v3.normalize.output_mod2.h5ad",
  min_counts_per_cell = 500,
  min_counts_per_gene = 500
)
## VIASH END

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

cat("Filtering genes\n")
sc$pp$filter_genes(ad_rna, min_counts = par$min_counts_per_gene)
sc$pp$filter_genes(ad_mod2, min_counts = par$min_counts_per_gene)

cat("Filtering cells\n")
mat <- cbind(ad_rna$X, ad_mod2$X)
filt <- rowSums(mat) > par$min_counts_per_cell
ad_rna <- ad_rna[filt, ]
ad_mod2 <- ad_mod2[filt, ]

cat("Storing output to '", par$output_rna, "'\n", sep = "")
zzz <- ad_rna$write_h5ad(par$output_rna, compression = "gzip")

cat("Storing output to '", par$output_mod2, "'\n", sep = "")
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")

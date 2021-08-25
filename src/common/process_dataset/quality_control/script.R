##
## TODO: Note that this script contains only very basic preprocessing.
##   This should be replaced at a later stage with something that makes more sense.
##

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(readr, warn.conflicts = FALSE, quietly = TRUE)

## VIASH START
par <- list(
  input_rna = "resources_test/common/test_resource.tmp.output_rna.h5ad",
  input_mod2 = "resources_test/common/test_resource.tmp.output_mod2.h5ad",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad",
  min_counts_per_cell = 10000,
  min_counts_per_gene = 150000,
  keep_genes = "src/common/resources/all_genes_tirosh.txt"
)
## VIASH END

cat("Reading h5ad file\n")
ad_rna <- anndata::read_h5ad(par$input_rna)
ad_mod2 <- anndata::read_h5ad(par$input_mod2)

obs <- ad_rna$obs
var_rna <- ad_rna$var
var_mod2 <- ad_mod2$var
mat_rna <- as(ad_rna$X, "CsparseMatrix")
mat_mod2 <- as(ad_mod2$X, "CsparseMatrix")

# read in which genes are supposed to be kept. multiple files are supported.
keep_genes <- unlist(lapply(par$keep_genes, readr::read_lines))

cat("Filtering genes and cells\n")
fil_rna_genes <- colSums(mat_rna) >= par$min_counts_per_gene | colnames(mat_rna) %in% keep_genes
fil_mod2_genes <- colSums(mat_mod2) >= par$min_counts_per_gene | colnames(mat_mod2) %in% keep_genes
fil_cells <- rowSums(mat_rna) >= par$min_counts_per_cell & rowSums(mat_mod2) >= par$min_counts_per_cell

obs <- obs[fil_cells, , drop = FALSE]
var_rna <- var_rna[fil_rna_genes, , drop = FALSE]
var_mod2 <- var_mod2[fil_mod2_genes, , drop = FALSE]
mat_rna <- as(mat_rna[fil_cells, fil_rna_genes, drop = FALSE], "CsparseMatrix")
mat_mod2 <- as(mat_mod2[fil_cells, fil_mod2_genes, drop = FALSE], "CsparseMatrix")

cat("Releveling factors\n")
if (ncol(obs) == 0) {
  obs <- NULL
} else {
  # relevel factors
  for (i in seq_len(ncol(obs))) {
    if (is.factor(obs[[i]])) {
      obs[[i]] <- forcats::fct_drop(obs[[i]])
    }
  }
  # relevel factors
  for (i in seq_len(ncol(var_rna))) {
    if (is.factor(var_rna[[i]])) {
      var_rna[[i]] <- forcats::fct_drop(var_rna[[i]])
    }
  }
  # relevel factors
  for (i in seq_len(ncol(var_mod2))) {
    if (is.factor(var_mod2[[i]])) {
      var_mod2[[i]] <- forcats::fct_drop(var_mod2[[i]])
    }
  }
}

cat("Writing mod1 data\n")
ad_rna_new <- anndata::AnnData(
  X = mat_rna,
  obs = obs,
  var = var_rna,
  uns = ad_rna$uns
)
print(ad_rna_new)
zzz <- ad_rna_new$write_h5ad(par$output_rna, compression = "gzip")

cat("Writing mod2 data\n")
ad_mod2_new <- anndata::AnnData(
  X = mat_mod2,
  obs = obs,
  var = var_mod2,
  uns = ad_mod2$uns
)
print(ad_mod2_new)
zzz <- ad_mod2_new$write_h5ad(par$output_mod2, compression = "gzip")

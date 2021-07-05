cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
sc <- reticulate::import("scanpy")

## VIASH START
par <- list(
  # input = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
  # output = "output1.h5ad"
  input = "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_raw_feature_bc_matrix.h5",
  output = "output2.h5ad"
)
## VIASH END

cat("Downloading file at '", par$input, "'\n", sep = "")
h5_tmp <- tempfile()
on.exit(file.remove(h5_tmp))
download.file(par$input, destfile = h5_tmp, quiet = TRUE)

cat("Reading 10x h5 file\n")
ad <- sc$read_10x_h5(h5_tmp, gex_only = FALSE)

cat("Setting dataset id\n")
ad$uns[["dataset_id"]] <- par$id

cat("Making var names unique\n")
ad$var_names_make_unique()

is_abseq <- ad$var$feature_types == "Antibody Capture"
if (any(is_abseq)) {
  cat("Storing Antibody Capture data in 'ad.obsm[\"protein\"]'\n")
  ab_counts <- ad$X[, is_abseq]
  ad <- ad[, !is_abseq]$copy()
  ad$obsm[["protein"]] <- as(ab_counts, "RsparseMatrix")
  ad$uns[["protein_varnames"]] <- colnames(ab_counts)
}

is_atacseq <- ad$var$feature_types == "Peaks"
if (any(is_atacseq)) {
  cat("Storing Peaks data in 'ad.obsm[\"chromatin\"]'\n")
  atac_counts <- ad$X[, is_atacseq]

  ad <- ad[, !is_atacseq]$copy()
  ad$obsm[["chromatin"]] <- as(atac_counts, "RsparseMatrix")
  ad$uns[["chromatin_varnames"]] <- colnames(atac_counts)
}


cat("Storing output to '", par$output, "'\n", sep = "")
ad$write_h5ad(par$output, compression = "gzip")

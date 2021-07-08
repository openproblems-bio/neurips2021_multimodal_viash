cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE, quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
sc <- reticulate::import("scanpy")

## VIASH START
par <- list(
  id = "lymph_node_lymphoma_14k",
  input = "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_raw_feature_bc_matrix.h5",
  output_rna = "output_rna.h5ad",
  output_mod2 = "output_mod2.h5ad"
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
is_atacseq <- ad$var$feature_types == "Peaks"

expect_true(
  any(is_abseq) != any(is_atacseq), 
  info = "Dataset should contain either antibody capture or ATAC peaks, not both"
)

if (any(is_abseq)) {
  cat("Processing Antibody data\n")
  ad_mod2 <- ad[, is_abseq]$copy()
  ad_mod2$uns[["modality"]] <- "Antibody"

  ad_mod1 <- ad[, !is_abseq]$copy()
  ad_mod1$uns[["modality"]] <- "RNA"
}

if (any(is_atacseq)) {
  cat("Processing ATAC data\n")
  ad_mod2 <- ad[, is_atacseq]$copy()
  ad_mod2$uns[["modality"]] <- "ATAC"

  ad_mod1 <- ad[, !is_atacseq]$copy()
  ad_mod1$uns[["modality"]] <- "RNA"
}


cat("Saving RNA data to '", par$output_rna, "'\n", sep = "")
zzz <- ad_mod1$write_h5ad(par$output_rna, compression = "gzip")

cat("Storing ", ad_mod2$uns[["modality"]], " data as '", par$output_mod2, "'\n", sep = "")
zzz <- ad_mod2$write_h5ad(par$output_mod2, compression = "gzip")
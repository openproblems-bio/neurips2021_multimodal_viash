library(testthat)
library(anndata)
#   # input = "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_raw_feature_bc_matrix.h5",
# input = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",


cat(">> Testing with citeseq dataset\n")
system(paste0(
  "./download_10x_dataset ",
  "--input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 ",
  "--output output1.h5ad ",
  "--id dataset1"
))

expect_true(file.exists("output1.h5ad"), info = "Output output1.h5ad was not created")

ad <- anndata::read_h5ad("output1.h5ad")

expect_is(ad$X, "sparseMatrix")
expect_gte(nrow(ad), 100)
expect_gte(ncol(ad), 100)

expect_is(ad$obsm[["protein"]], "sparseMatrix")
expect_gte(nrow(ad$obsm[["protein"]]), 100)
expect_gte(ncol(ad$obsm[["protein"]]), 10)
expect_equal(length(ad$uns[["protein_varnames"]]), ncol(ad$obsm[["protein"]]))

cat(">> Testing with atacseq dataset\n")
system(paste0(
  "./download_10x_dataset ",
  "--input https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_raw_feature_bc_matrix.h5 ",
  "--output output2.h5ad ",
  "--id dataset2"
))

expect_true(file.exists("output2.h5ad"), info = "Output output2.h5ad was not created")

ad <- anndata::read_h5ad("output2.h5ad")
expect_is(ad$X, "sparseMatrix")
expect_gte(nrow(ad), 100)
expect_gte(ncol(ad), 100)

expect_is(ad$obsm[["chromatin"]], "sparseMatrix")
expect_gte(nrow(ad$obsm[["chromatin"]]), 100)
expect_gte(ncol(ad$obsm[["chromatin"]]), 10)
expect_equal(length(ad$uns[["chromatin_varnames"]]), ncol(ad$obsm[["chromatin"]]))



cat(">> All tests passed successfully!\n")

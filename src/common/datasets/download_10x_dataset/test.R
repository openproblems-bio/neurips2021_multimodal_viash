library(testthat)
library(anndata)

cat(">> Testing with citeseq dataset\n")
system(paste0(
  "./download_10x_dataset ",
  "--input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 ",
  "--organism human ",
  "--output_rna output1_rna.h5ad ",
  "--output_mod2 output1_mod2.h5ad ",
  "--id dataset1"
))

cat(">> Checking whether output files exist\n")
expect_true(file.exists("output1_rna.h5ad"), info = "Output output1_rna.h5ad was not created")
expect_true(file.exists("output1_mod2.h5ad"), info = "Output output1_mod2.h5ad was not created")

cat(">> Reading output files\n")
ad1 <- anndata::read_h5ad("output1_rna.h5ad")
ad2 <- anndata::read_h5ad("output1_mod2.h5ad")

cat(">> Checking RNA data\n")
expect_is(ad1$X, "sparseMatrix")
expect_gte(nrow(ad1), 100)
expect_gte(ncol(ad1), 100)
expect_equal(ad1$uns[["dataset_id"]], "dataset1")

cat(">> Checking antibody data\n")
expect_is(ad2$X, "sparseMatrix")
expect_gte(nrow(ad2), 100)
expect_gte(ncol(ad2), 10)
expect_equal(ad2$uns[["dataset_id"]], "dataset1")

expect_equal(rownames(ad1), rownames(ad2))

cat(">> Testing with atacseq dataset\n")
system(paste0(
  "./download_10x_dataset ",
  "--input https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_raw_feature_bc_matrix.h5 ",
  "--output_rna output2_rna.h5ad ",
  "--output_mod2 output2_mod2.h5ad ",
  "--id dataset2"
))

cat(">> Checking whether output files exist\n")
expect_true(file.exists("output2_rna.h5ad"), info = "Output output2_rna.h5ad was not created")
expect_true(file.exists("output2_mod2.h5ad"), info = "Output output2_mod2.h5ad was not created")

cat(">> Reading output files\n")
ad1 <- anndata::read_h5ad("output2_rna.h5ad")
ad2 <- anndata::read_h5ad("output2_mod2.h5ad")

cat(">> Checking RNA data\n")
expect_is(ad1$X, "sparseMatrix")
expect_gte(nrow(ad1), 100)
expect_gte(ncol(ad1), 100)
expect_equal(ad1$uns[["dataset_id"]], "dataset2")
expect_equal(ad1$uns[["organism"]], "human")

cat(">> Checking antibody data\n")
expect_is(ad2$X, "sparseMatrix")
expect_gte(nrow(ad2), 100)
expect_gte(ncol(ad2), 100)
expect_equal(ad2$uns[["dataset_id"]], "dataset2")

expect_equal(rownames(ad1), rownames(ad2))


cat(">> All tests passed successfully!\n")

library(testthat, quietly = TRUE, warn.conflicts = FALSE)
library(anndata, quietly = TRUE, warn.conflicts = FALSE)

cat(">> Run normalize component\n")
system(paste0(
  "./normalize ",
  "--input_rna resources_test/common/pbmc_1k_protein_v3.normalize.output_rna.h5ad ",
  "--input_mod2 resources_test/common/pbmc_1k_protein_v3.normalize.output_mod2.h5ad ",
  "--output_rna output_rna.h5ad ",
  "--output_mod2 output_mod2.h5ad ",
  "--min_counts_per_gene 1200 ",
  "--min_counts_per_cell 1200"
))

cat(">> Checking whether output files were created\n")
expect_true(file.exists("output_rna.h5ad"), info = "Output output_rna.h5ad was not created")
expect_true(file.exists("output_mod2.h5ad"), info = "Output output_mod2.h5ad was not created")

cat(">> Reading in h5ad files\n")
ad_in1 <- anndata::read_h5ad("resources_test/common/pbmc_1k_protein_v3.normalize.output_rna.h5ad")
ad_in2 <- anndata::read_h5ad("resources_test/common/pbmc_1k_protein_v3.normalize.output_mod2.h5ad")
ad_out1 <- anndata::read_h5ad("output_rna.h5ad")
ad_out2 <- anndata::read_h5ad("output_mod2.h5ad")

cat(">> Checking rna output\n")
expect_is(ad_out1$X, "sparseMatrix")
expect_lte(nrow(ad_out1), nrow(ad_in1))
expect_lte(ncol(ad_out1), ncol(ad_in1))
expect_equal(ad_out1$uns[["dataset_id"]], ad_in1$uns[["dataset_id"]])
expect_equal(unique(ad_out1$var[["feature_types"]]), unique(ad_in1$var[["feature_types"]]))

cat(">> Checking mod2 output\n")
expect_is(ad_out2$X, "sparseMatrix")
expect_lte(nrow(ad_out2), nrow(ad_in2))
expect_lte(ncol(ad_out2), ncol(ad_in2))
expect_equal(ad_out2$uns[["dataset_id"]], ad_in2$uns[["dataset_id"]])
expect_equal(unique(ad_out2$var[["feature_types"]]), unique(ad_in2$var[["feature_types"]]))

cat(">> All tests passed successfully!\n")

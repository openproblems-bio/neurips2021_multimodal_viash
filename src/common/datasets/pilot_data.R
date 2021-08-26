library(tidyverse)
library(Matrix)
library(anndata)

# sync to local folder
system("aws s3 sync --profile op2 s3://openproblems-bio/neurips2021/processed/atac/ output/private_datasets/common_unproc/atac/")

dir.create("output/private_datasets/common/openproblems_bmmc_multiome", recursive = TRUE)
# process private datasets
ad1 <- read_h5ad("output/private_datasets/common_unproc/atac/Multiome_gex_processed.h5ad")
ad1$uns[["dataset_id"]] <- "openproblems_bmmc_multiome"
ad1$obs[["is_train"]] <- ad1$obs[["batch"]] != "s1d2"
ad1$X <- as(ad1$X, "CsparseMatrix")
ad1$write_h5ad("output/private_datasets/common/openproblems_bmmc_multiome/openproblems_bmmc_multiome.manual_formatting.output_rna.h5ad", compression = "gzip")

ad2 <- read_h5ad("output/private_datasets/common_unproc/atac/Multiome_atac_processed_gact.h5ad")
ad2$uns[["dataset_id"]] <- "openproblems_bmmc_multiome"
ad2$obs[["is_train"]] <- ad2$obs[["batch"]] != "s1d2"
ad2$X <- as(ad2$X, "CsparseMatrix")
ad2$write_h5ad("output/private_datasets/common/openproblems_bmmc_multiome/openproblems_bmmc_multiome.manual_formatting.output_mod2.h5ad", compression = "gzip")

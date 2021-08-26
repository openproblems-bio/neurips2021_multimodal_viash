library(tidyverse)
library(Matrix)
library(anndata)

# sync to local folder
# system("aws s3 sync --profile op2 s3://openproblems-bio/neurips2021/processed/atac/ output/manual_formatting/atac/")

# switch over to public at some point
output_dir <- "output/public_datasets/common/openproblems_bmmc_multiome"

dir.create(output_dir, recursive = TRUE)

# read data
ad1 <- read_h5ad("output/manual_formatting/atac/Multiome_gex_processed.h5ad")
ad2 <- read_h5ad("output/manual_formatting/atac/Multiome_atac_processed.h5ad")

# process gex
ad1$uns[["dataset_id"]] <- "openproblems_bmmc_multiome"
ad1$uns[["organism"]] <- "human"
ad1$obs[["is_train"]] <- ad1$obs[["batch"]] != "s1d2"
ad1$X <- as(ad1$layers[["counts"]], "CsparseMatrix")
ad1$layers <- NULL

# process atac
ad2$uns[["dataset_id"]] <- "openproblems_bmmc_multiome"
ad2$uns[["organism"]] <- "human"
ad2$obs[["is_train"]] <- ad2$obs[["batch"]] != "s1d2"
ad2$X@x <- (ad2$X@x > 0) + 0

# copy over data
ad1$obs[["pseudotime_order_ATAC"]] <- ad2$obs[["pseudotime_order_ATAC"]]
ad2$obs[["pseudotime_order_GEX"]] <- ad1$obs[["pseudotime_order_GEX"]]

# write to file
ad1$write_h5ad(paste0(output_dir, "/openproblems_bmmc_multiome.manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2$write_h5ad(paste0(output_dir, "/openproblems_bmmc_multiome.manual_formatting.output_mod2.h5ad"), compression = "gzip")

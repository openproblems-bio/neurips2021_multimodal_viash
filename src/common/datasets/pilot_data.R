library(tidyverse)
library(Matrix)
library(anndata)

# sync to local folder
# system("aws s3 sync s3://openproblems-bio/public/ output/manual_formatting/ --no-sign-request")

#############################
# Process multiome

output_dir <- "output/public_datasets/common/openproblems_bmmc_multiome"

dir.create(output_dir, recursive = TRUE)

# read data
ad1 <- read_h5ad("output/manual_formatting/multiome/Multiome_GEX_processed.h5ad")
ad2 <- read_h5ad("output/manual_formatting/multiome/Multiome_ATAC_processed.h5ad")

# check reads
ad1$layers[["counts"]][1:10, 1:10]
ad2$X[1:10, 1:10]
head(ad1$obs)
head(ad2$obs)
head(ad1$var)
head(ad2$var)

# process gex
# ad1$uns[["dataset_id"]] <- "openproblems_bmmc_multiome"
# ad1$uns[["organism"]] <- "human"
ad1$uns[["batch_type"]] <- "real"
ad1$obs[["is_train"]] <- ad1$obs[["batch"]] != "s1d2"
ad1$X <- as(ad1$layers[["counts"]], "CsparseMatrix")
ad1$layers <- NULL

# process atac
# ad2$uns[["dataset_id"]] <- "openproblems_bmmc_multiome"
# ad2$uns[["organism"]] <- "human"
ad2$uns[["batch_type"]] <- "real"
ad2$obs[["is_train"]] <- ad2$obs[["batch"]] != "s1d2"
ad2$X@x <- (ad2$X@x > 0) + 0

# copy over data
ad1$obs[["pseudotime_order_ATAC"]] <- ad2$obs[["pseudotime_order_ATAC"]]
ad2$obs[["pseudotime_order_GEX"]] <- ad1$obs[["pseudotime_order_GEX"]]

# write to file
ad1$write_h5ad(paste0(output_dir, "/", ad1$uns[["dataset_id"]], ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2$write_h5ad(paste0(output_dir, "/", ad1$uns[["dataset_id"]], ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


#############################
# Process cite

output_dir <- "output/public_datasets/common/openproblems_bmmc_cite"

dir.create(output_dir, recursive = TRUE)

# read data
ad1 <- read_h5ad("output/manual_formatting/cite/CITE_GEX_processed.h5ad")
ad2 <- read_h5ad("output/manual_formatting/cite/CITE_ADT_processed.h5ad")

# check reads
ad1$layers[["counts"]][1:10, 1:10]
ad2$X[1:10, 1:10]
head(ad1$obs)
head(ad2$obs)
head(ad1$var)
head(ad2$var)

# process gex
ad1$uns[["batch_type"]] <- "real"
ad1$obs[["is_train"]] <- ad1$obs[["batch"]] != "s1d2"
ad1$X <- as(ad1$layers[["counts"]], "CsparseMatrix")
ad1$layers <- NULL

# process atac
ad2$uns[["batch_type"]] <- "real"
ad2$obs[["is_train"]] <- ad2$obs[["batch"]] != "s1d2"

# copy over data
ad1$obs[["pseudotime_order_ADT"]] <- ad2$obs[["pseudotime_order_ADT"]]
ad2$obs[["pseudotime_order_GEX"]] <- ad1$obs[["pseudotime_order_GEX"]]

# write to file
ad1$write_h5ad(paste0(output_dir, "/", ad1$uns[["dataset_id"]], ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2$write_h5ad(paste0(output_dir, "/", ad1$uns[["dataset_id"]], ".manual_formatting.output_mod2.h5ad"), compression = "gzip")

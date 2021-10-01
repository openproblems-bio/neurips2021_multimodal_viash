library(tidyverse)
library(Matrix)
library(anndata)
library(assertthat)

set.seed(1)


resources_test <- "resources_test/common/"
output_dir <- "output/datasets_loocv/common/"


#############################
# Process multiome

# read data
ad1 <- read_h5ad("output/manual_formatting/multiome/Multiome_gex_processed_train.h5ad")
ad2 <- read_h5ad("output/manual_formatting/multiome/Multiome_atac_processed_train.h5ad")
adid <- "openproblems_bmmc_multiome"
assert_that(all(rownames(ad1) == rownames(ad2)))

levels(ad1$obs$batch)
assert_that(length(setdiff(ad1$obs$batch, c(train, valid, backup_test))) == 0)
assert_that(!any(backup_test %in% levels(ad1$obs$batch)))
assert_that(!any(backup_test %in% levels(ad2$obs$batch)))
sort(setdiff(c(train, valid, backup_test), unique(ad1$obs$batch)))

# check reads
ad1$layers[["counts"]][1:10, 1:10]
ad2$X[1:10, 1:10]
head(ad1$obs)
head(ad2$obs)
head(ad1$var)
head(ad2$var)

# process gex
ad1$X <- as(ad1$layers[["log_norm"]], "CsparseMatrix")
ad1$layers <- list(counts = as(ad1$layers[["counts"]], "CsparseMatrix"))

# process atac
ad2$X <- as(ad2$X, "CsparseMatrix")
ad2$layers <- list(counts = as(ad2$layers[["counts"]], "CsparseMatrix"))

assert_that(max(ad1$X) < 100)
assert_that(max(ad2$X) < 100)

# copy over data
ad1$obs[["pseudotime_order_ATAC"]] <- ad2$obs[["pseudotime_order_ATAC"]]
ad2$obs[["pseudotime_order_GEX"]] <- ad1$obs[["pseudotime_order_GEX"]]

#############
# create datasets
abatches <- levels(ad1$obs$batch)

for (loo_batch in abatches) {
  ad1_loocv <- ad1$copy()
  ad2_loocv <- ad2$copy()

  ad1_loocv$obs[["is_train"]] <- ad1_loocv$obs[["batch"]] != loo_batch
  ad2_loocv$obs[["is_train"]] <- ad1_loocv$obs[["batch"]] != loo_batch
  
  adid_loocv <- paste0(adid, "_loocv_", loo_batch)
  ad1_loocv$uns[["dataset_id"]] <- adid_loocv
  ad2_loocv$uns[["dataset_id"]] <- adid_loocv
  
  dir.create(paste0(output_dir, adid_loocv), recursive = TRUE)
  ad1_loocv$write_h5ad(paste0(output_dir, adid_loocv, "/", adid_loocv, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
  ad2_loocv$write_h5ad(paste0(output_dir, adid_loocv, "/", adid_loocv, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")
}


#############################
# Process cite

# read data
bd1 <- read_h5ad("output/manual_formatting/cite/Cite_gex_processed_train.h5ad")
bd2 <- read_h5ad("output/manual_formatting/cite/Cite_adt_processed_train.h5ad")
bdid <- "openproblems_bmmc_cite"
assert_that(all(rownames(bd1) == rownames(bd2)))

levels(bd1$obs$batch)
assert_that(!any(backup_test %in% levels(bd1$obs$batch)))
assert_that(!any(backup_test %in% levels(bd2$obs$batch)))
testthat::expect_length(setdiff(bd1$obs$batch, c(train, valid, backup_test)), 0)
setdiff(c(train, valid, backup_test), unique(bd1$obs$batch))

# check reads
bd1$layers[["counts"]][1:10, 1:10]
bd2$layers[["counts"]][1:10, 1:10]
head(bd1$obs)
head(bd2$obs)
head(bd1$var)
head(bd2$var)

# process gex
bd1$X <- as(bd1$layers[["log_norm"]], "CsparseMatrix")
bd1$layers <- list(counts = as(bd1$layers[["counts"]], "CsparseMatrix"))

# process atac
bd2$X <- as(bd2$X, "CsparseMatrix")
bd2$layers <- list(counts = as(bd2$layers[["counts"]], "CsparseMatrix"))

assert_that(max(bd1$X) < 100)
assert_that(max(bd2$X) < 100)

# copy over data
bd1$obs[["pseudotime_order_ADT"]] <- bd2$obs[["pseudotime_order_ADT"]]
bd2$obs[["pseudotime_order_GEX"]] <- bd1$obs[["pseudotime_order_GEX"]]

#############
# create loocv datasets

bbatches <- levels(bd1$obs$batch)

for (loo_batch in bbatches) {
  bd1_loocv <- bd1$copy()
  bd2_loocv <- bd2$copy()
  
  bd1_loocv$obs[["is_train"]] <- bd1_loocv$obs[["batch"]] != loo_batch
  bd2_loocv$obs[["is_train"]] <- bd1_loocv$obs[["batch"]] != loo_batch
  
  bdid_loocv <- paste0(bdid, "_loocv_", loo_batch)
  bd1_loocv$uns[["dataset_id"]] <- bdid_loocv
  bd2_loocv$uns[["dataset_id"]] <- bdid_loocv
  
  dir.create(paste0(output_dir, bdid_loocv), recursive = TRUE)
  bd1_loocv$write_h5ad(paste0(output_dir, bdid_loocv, "/", bdid_loocv, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
  bd2_loocv$write_h5ad(paste0(output_dir, bdid_loocv, "/", bdid_loocv, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")
}
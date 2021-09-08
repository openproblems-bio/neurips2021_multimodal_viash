library(tidyverse)
library(Matrix)
library(anndata)

set.seed(1)

# sync to local folder
# system("aws s3 sync s3://openproblems-bio/public/explore/ output/manual_formatting/ --profile op2")

starter_dev <- c("s1d1", "s1d2")
train <- c("s1d1", "s2d1", "s2d4", "s3d6", "s3d1")
valid <- c("s1d2", "s3d7")
backup_test <- c("s1d3", "s2d5")

resources_test <- "resources_test/common/"
output_dir <- "output/datasets/common/"

#############################
# Process multiome

# read data
ad1 <- read_h5ad("output/manual_formatting/multiome/Multiome_GEX_processed.training.h5ad")
ad2 <- read_h5ad("output/manual_formatting/multiome/Multiome_ATAC_processed.training.h5ad")
adid <- ad1$uns[["dataset_id"]]

unique(ad1$obs$batch)
testthat::expect_length(setdiff(ad1$obs$batch, c(train, valid, backup_test)), 0)
setdiff(c(train, valid, backup_test), unique(ad1$obs$batch))

# check reads
ad1$layers[["counts"]][1:10, 1:10]
ad2$X[1:10, 1:10]
head(ad1$obs)
head(ad2$obs)
head(ad1$var)
head(ad2$var)

# process gex
ad1$uns[["batch_type"]] <- "real"
ad1$X <- as(ad1$layers[["counts"]], "CsparseMatrix")
ad1$layers <- NULL

# process atac
ad2$uns[["batch_type"]] <- "real"
ad2$X@x <- (ad2$X@x > 0) + 0

# copy over data
ad1$obs[["pseudotime_order_ATAC"]] <- ad2$obs[["pseudotime_order_ATAC"]]
ad2$obs[["pseudotime_order_GEX"]] <- ad1$obs[["pseudotime_order_GEX"]]

#############
# create phase 1

ad1_phase1 <- ad1[ad1$obs$batch %in% c(train, valid), ]$copy()
ad2_phase1 <- ad2[ad2$obs$batch %in% c(train, valid), ]$copy()

ad1_phase1$obs[["is_train"]] <- ad1_phase1$obs[["batch"]] %in% train
ad2_phase1$obs[["is_train"]] <- ad1_phase1$obs[["batch"]] %in% train

adid_phase1 <- paste0(adid, "_phase1")
ad1_phase1$uns[["dataset_id"]] <- adid_phase1
ad2_phase1$uns[["dataset_id"]] <- adid_phase1

dir.create(paste0(output_dir, adid_phase1), recursive = TRUE)
ad1_phase1$write_h5ad(paste0(output_dir, adid_phase1, "/", adid_phase1, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_phase1$write_h5ad(paste0(output_dir, adid_phase1, "/", adid_phase1, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


#############
# create phase 2

ad1_phase2 <- ad1$copy()
ad2_phase2 <- ad2$copy()

ad1_phase2$obs[["is_train"]] <- ad1_phase2$obs[["batch"]] %in% c(train, valid)
ad2_phase2$obs[["is_train"]] <- ad1_phase2$obs[["batch"]] %in% c(train, valid)

adid_phase2 <- paste0(adid, "_phase2")
ad1_phase2$uns[["dataset_id"]] <- adid_phase2
ad2_phase2$uns[["dataset_id"]] <- adid_phase2

dir.create(paste0(output_dir, adid_phase2), recursive = TRUE)
ad1_phase2$write_h5ad(paste0(output_dir, adid_phase2, "/", adid_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_phase2$write_h5ad(paste0(output_dir, adid_phase2, "/", adid_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


#############
# create openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter

starter_cell <- which(ad1$obs$batch %in% starter_dev)
row_sel <- sort(sample(starter_cell, 500))
col_sel1 <- seq_len(ncol(ad1)) %in% sample.int(ncol(ad1), 600, prob = Matrix::colSums(ad1$X > 0))
col_sel2 <- seq_len(ncol(ad2)) %in% sample.int(ncol(ad2), 600, prob = Matrix::colSums(ad2$X > 0))

ad1_starter <- ad1[row_sel, col_sel1]$copy()
ad2_starter <- ad2[row_sel, col_sel2]$copy()

ad1_starter$obsm <- NULL
ad2_starter$obsm <- NULL

ad1_starter$obs[["is_train"]] <- ad1_starter$obs[["batch"]] %in% train
ad2_starter$obs[["is_train"]] <- ad2_starter$obs[["batch"]] %in% train

adid_starter <- paste0(adid, "_starter")
ad1_starter$uns[["dataset_id"]] <- adid_starter
ad2_starter$uns[["dataset_id"]] <- adid_starter

dir.create(paste0(resources_test, adid_starter), recursive = TRUE)
ad1_starter$write_h5ad(paste0(resources_test, adid_starter, "/", adid_starter, ".output_rna.h5ad"), compression = "gzip")
ad2_starter$write_h5ad(paste0(resources_test, adid_starter, "/", adid_starter, ".output_mod2.h5ad"), compression = "gzip")





#############################
# Process cite

# read data
bd1 <- read_h5ad("output/manual_formatting/cite/CITE_GEX_processed.training.h5ad")
bd2 <- read_h5ad("output/manual_formatting/cite/CITE_ADT_processed.training.h5ad")
bdid <- bd1$uns[["dataset_id"]]

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
bd1$uns[["batch_type"]] <- "real"
bd1$X <- as(bd1$layers[["counts"]], "CsparseMatrix")
bd1$layers <- NULL

# process atac
bd2$uns[["batch_type"]] <- "real"
bd2$X <- as(bd2$layers[["counts"]], "CsparseMatrix")
bd2$layers <- NULL

# copy over data
bd1$obs[["pseudotime_order_ADT"]] <- bd2$obs[["pseudotime_order_ADT"]]
bd2$obs[["pseudotime_order_GEX"]] <- bd1$obs[["pseudotime_order_GEX"]]

#############
# create phase 1

bd1_phase1 <- bd1[bd1$obs$batch %in% c(train, valid), ]$copy()
bd2_phase1 <- bd2[bd2$obs$batch %in% c(train, valid), ]$copy()

bd1_phase1$obs[["is_train"]] <- bd1_phase1$obs[["batch"]] %in% train
bd2_phase1$obs[["is_train"]] <- bd1_phase1$obs[["batch"]] %in% train

bdid_phase1 <- paste0(bdid, "_phase1")
bd1_phase1$uns[["dataset_id"]] <- bdid_phase1
bd2_phase1$uns[["dataset_id"]] <- bdid_phase1

dir.create(paste0(output_dir, bdid_phase1), recursive = TRUE)
bd1_phase1$write_h5ad(paste0(output_dir, bdid_phase1, "/", bdid_phase1, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
bd2_phase1$write_h5ad(paste0(output_dir, bdid_phase1, "/", bdid_phase1, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


# #############
# # create phase 2

# bd1_phase2 <- bd1$copy()
# bd2_phase2 <- bd2$copy()

# bd1_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)
# bd2_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)

# bdid_phase2 <- paste0(bdid, "_phase2")
# bd1_phase2$uns[["dataset_id"]] <- bdid_phase2
# bd2_phase2$uns[["dataset_id"]] <- bdid_phase2

# dir.create(paste0(output_dir, bdid_phase2), recursive = TRUE)
# bd1_phase2$write_h5ad(paste0(output_dir, bdid_phase2, "/", bdid_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
# bd2_phase2$write_h5ad(paste0(output_dir, bdid_phase2, "/", bdid_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")



#############
# create openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter

starter_cell <- which(bd1$obs$batch %in% starter_dev)
row_sel <- sort(sample(starter_cell, 600))
col_sel1 <- seq_len(ncol(bd1)) %in% sample.int(ncol(bd1), 1000, prob = Matrix::colSums(bd1$X > 0))

bd1_starter <- bd1[row_sel, col_sel1]$copy()
bd2_starter <- bd2[row_sel, ]$copy()

bd1_starter$obs[["is_train"]] <- bd1_starter$obs[["batch"]] %in% train
bd2_starter$obs[["is_train"]] <- bd2_starter$obs[["batch"]] %in% train

bdid_starter <- paste0(bdid, "_starter")
bd1_starter$uns[["dataset_id"]] <- bdid_starter
bd2_starter$uns[["dataset_id"]] <- bdid_starter

bd1_starter$obsm <- NULL
bd2_starter$obsm <- NULL

dir.create(paste0(resources_test, bdid_starter), recursive = TRUE)
bd1_starter$write_h5ad(paste0(resources_test, bdid_starter, "/", bdid_starter, ".output_rna.h5ad"), compression = "gzip")
bd2_starter$write_h5ad(paste0(resources_test, bdid_starter, "/", bdid_starter, ".output_mod2.h5ad"), compression = "gzip")

library(tidyverse)
library(Matrix)
library(anndata)

set.seed(1)

# sync to local folder
# system("aws s3 sync s3://openproblems-bio/public/ output/manual_formatting/ --profile op2")

train <- c("s1d1", "s1d2", "s2d4", "s3d1")
valid <- c("s2d1", "s3d6")
backup_test <- c("s1d3", "s2d5", "s3d7")
output_dir <- "output/datasets/common/"

#############################
# Process multiome

# read data
ad1 <- read_h5ad("output/manual_formatting/multiome/Multiome_GEX_processed.h5ad")
ad2 <- read_h5ad("output/manual_formatting/multiome/Multiome_ATAC_processed.h5ad")
dataset_id <- ad1$uns[["dataset_id"]]

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

iid_sel <- seq_len(nrow(ad1)) %in% sample.int(nrow(ad1), 0.1 * nrow(ad1))

table(iid_sel)
table(iid_sel, ad1$obs$batch)

#############
# create phase 1

ad1_phase1 <- ad1[!iid_sel & ad1$obs$batch %in% c(train, valid), ]$copy()
ad2_phase1 <- ad2[!iid_sel & ad2$obs$batch %in% c(train, valid), ]$copy()

ad1_phase1$obs[["is_train"]] <- ad1_phase1$obs[["batch"]] %in% train
ad2_phase1$obs[["is_train"]] <- ad1_phase1$obs[["batch"]] %in% train

did_phase1 <- paste0(dataset_id, "_phase1")
ad1_phase1$uns[["dataset_id"]] <- did_phase1
ad2_phase1$uns[["dataset_id"]] <- did_phase1

dir.create(paste0(output_dir, did_phase1), recursive = TRUE)
ad1_phase1$write_h5ad(paste0(output_dir, did_phase1, "/", did_phase1, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_phase1$write_h5ad(paste0(output_dir, did_phase1, "/", did_phase1, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


#############
# create phase 2

ad1_phase2 <- ad1$copy()
ad2_phase2 <- ad2$copy()

ad1_phase2$obs[["is_train"]] <- ad1_phase2$obs[["batch"]] %in% c(train, valid)
ad2_phase2$obs[["is_train"]] <- ad1_phase2$obs[["batch"]] %in% c(train, valid)

did_phase2 <- paste0(dataset_id, "_phase2")
ad1_phase2$uns[["dataset_id"]] <- did_phase2
ad2_phase2$uns[["dataset_id"]] <- did_phase2

dir.create(paste0(output_dir, did_phase2), recursive = TRUE)
ad1_phase2$write_h5ad(paste0(output_dir, did_phase2, "/", did_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_phase2$write_h5ad(paste0(output_dir, did_phase2, "/", did_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


#############
# create iid

ad1_iid <- ad1$copy()
ad2_iid <- ad2$copy()

ad1_iid$obs[["is_train"]] <- !iid_sel
ad2_iid$obs[["is_train"]] <- !iid_sel

did_iid <- paste0(dataset_id, "_iid")
ad1_iid$uns[["dataset_id"]] <- did_iid
ad2_iid$uns[["dataset_id"]] <- did_iid

dir.create(paste0(output_dir, did_iid), recursive = TRUE)
ad1_iid$write_h5ad(paste0(output_dir, did_iid, "/", did_iid, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_iid$write_h5ad(paste0(output_dir, did_iid, "/", did_iid, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")





#############################
# Process cite

# temporary setup until more samples arrive
train <- c("s1d1")
valid <- c("s1d2")
backup_test <- NULL


# read data
bd1 <- read_h5ad("output/manual_formatting/cite/CITE_GEX_processed.h5ad")
bd2 <- read_h5ad("output/manual_formatting/cite/CITE_ADT_processed.h5ad")
dataset_id <- bd1$uns[["dataset_id"]]

testthat::expect_length(setdiff(bd1$obs$batch, c(train, valid, backup_test)), 0)
setdiff(c(train, valid, backup_test), unique(bd1$obs$batch))

# check reads
bd1$layers[["counts"]][1:10, 1:10]
bd2$X[1:10, 1:10]
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

# copy over data
bd1$obs[["pseudotime_order_ADT"]] <- bd2$obs[["pseudotime_order_ADT"]]
bd2$obs[["pseudotime_order_GEX"]] <- bd1$obs[["pseudotime_order_GEX"]]

iid_sel <- seq_len(nrow(bd1)) %in% sample.int(nrow(bd1), 0.1 * nrow(bd1))

#############
# create phase 1

bd1_phase1 <- bd1[!iid_sel & bd1$obs$batch %in% c(train, valid), ]$copy()
bd2_phase1 <- bd2[!iid_sel & bd2$obs$batch %in% c(train, valid), ]$copy()

bd1_phase1$obs[["is_train"]] <- bd1_phase1$obs[["batch"]] %in% train
bd2_phase1$obs[["is_train"]] <- bd1_phase1$obs[["batch"]] %in% train

did_phase1 <- paste0(dataset_id, "_phase1")
bd1_phase1$uns[["dataset_id"]] <- did_phase1
bd2_phase1$uns[["dataset_id"]] <- did_phase1

dir.create(paste0(output_dir, did_phase1), recursive = TRUE)
bd1_phase1$write_h5ad(paste0(output_dir, did_phase1, "/", did_phase1, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
bd2_phase1$write_h5ad(paste0(output_dir, did_phase1, "/", did_phase1, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


# #############
# # create phase 2

# bd1_phase2 <- bd1$copy()
# bd2_phase2 <- bd2$copy()

# bd1_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)
# bd2_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)

# did_phase2 <- paste0(dataset_id, "_phase2")
# bd1_phase2$uns[["dataset_id"]] <- did_phase2
# bd2_phase2$uns[["dataset_id"]] <- did_phase2

# dir.create(paste0(output_dir, did_phase2), recursive = TRUE)
# bd1_phase2$write_h5ad(paste0(output_dir, did_phase2, "/", did_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
# bd2_phase2$write_h5ad(paste0(output_dir, did_phase2, "/", did_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


#############
# create iid

bd1_iid <- bd1$copy()
bd2_iid <- bd2$copy()

bd1_iid$obs[["is_train"]] <- !iid_sel
bd2_iid$obs[["is_train"]] <- !iid_sel

did_iid <- paste0(dataset_id, "_iid")
bd1_iid$uns[["dataset_id"]] <- did_iid
bd2_iid$uns[["dataset_id"]] <- did_iid

dir.create(paste0(output_dir, did_iid), recursive = TRUE)
bd1_iid$write_h5ad(paste0(output_dir, did_iid, "/", did_iid, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
bd2_iid$write_h5ad(paste0(output_dir, did_iid, "/", did_iid, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")


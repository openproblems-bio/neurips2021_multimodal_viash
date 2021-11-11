library(tidyverse)
library(Matrix)
library(anndata)
library(assertthat)

set.seed(1)

# sync to local folder, run in bash
# aws s3 sync s3://openproblems-bio/neurips2021/processed/phase2/ output/manual_formatting_2021-11-08/ --profile op2 --delete
# aws s3 sync s3://openproblems-bio/public/phase1-data/ output/datasets/ --profile op2 --delete --dryrun

starter_dev <- c("s1d1", "s1d2")
train <- c("s1d1", "s1d3", "s2d1", "s2d4", "s2d5", "s3d3", "s3d6", "s3d10")
valid <- c("s1d2", "s3d7")
test <- c("s4d1", "s4d8", "s4d9")

resources_test <- "resources_test/common/"
output_dir <- "output/datasets_2021-11-08/common/"


#############################
# Process multiome

# read data
ad1tr <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_gex_processed_train.h5ad")
ad2tr <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_atac_processed_train.h5ad")
ad1te <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_gex_processed_test_donors.h5ad")
ad2te <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_atac_processed_test_donors.h5ad")
ad1 <- anndata::concat(list(ad1tr, ad1te))
ad1$obs <- ad1$obs %>% mutate_if(is.character, factor)
ad1$var <- ad1tr$var %>% select(-contains("-s4d"))
ad1$uns <- ad1tr$uns
ad2 <- anndata::concat(list(ad2tr, ad2te))
ad2$obs <- ad2$obs %>% mutate_if(is.character, factor)
ad2$uns <- ad2tr$uns
ad2$var <- ad2tr$var %>% select(-contains("-s4d"))
adid <- "openproblems_bmmc_multiome"
assert_that(all(rownames(ad1) == rownames(ad2)))

unique(ad1$obs$batch)
assert_that(length(setdiff(ad1$obs$batch, c(train, valid, test))) == 0)
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
bd1tr <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_gex_processed_train.h5ad")
bd2tr <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_adt_processed_train.h5ad")
bd1te <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_gex_processed_test_donors.h5ad")
bd2te <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_adt_processed_test_donors.h5ad")
bd1 <- anndata::concat(list(bd1tr, bd1te))
bd1$obs <- bd1$obs %>% mutate_if(is.character, factor)
bd1$var <- bd1tr$var %>% select(-contains("-s4d"))
bd1$uns <- bd1tr$uns
bd2 <- anndata::concat(list(bd2tr, bd2te))
bd2$obs <- bd2$obs %>% mutate_if(is.character, factor)
bd2$uns <- bd2tr$uns
bd2$var <- bd2tr$var %>% select(-contains("-s4d"))
bdid <- "openproblems_bmmc_cite"
assert_that(all(rownames(bd1) == rownames(bd2)))

levels(bd1$obs$batch)
testthat::expect_length(setdiff(bd1$obs$batch, c(train, valid, test)), 0)
setdiff(c(train, valid, test), unique(bd1$obs$batch))

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


#############
# create phase 2

bd1_phase2 <- bd1$copy()
bd2_phase2 <- bd2$copy()

bd1_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)
bd2_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)

bdid_phase2 <- paste0(bdid, "_phase2")
bd1_phase2$uns[["dataset_id"]] <- bdid_phase2
bd2_phase2$uns[["dataset_id"]] <- bdid_phase2

dir.create(paste0(output_dir, bdid_phase2), recursive = TRUE)
bd1_phase2$write_h5ad(paste0(output_dir, bdid_phase2, "/", bdid_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
bd2_phase2$write_h5ad(paste0(output_dir, bdid_phase2, "/", bdid_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")



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

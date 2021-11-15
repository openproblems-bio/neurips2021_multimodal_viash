library(tidyverse)
library(Matrix)
library(anndata)
library(assertthat)

set.seed(1)

# sync to local folder, run in bash
# aws s3 sync s3://openproblems-bio/neurips2021/processed/phase2/ output/manual_formatting_2021-11-08/ --profile op2 --delete --dryrun
# aws s3 sync s3://openproblems-bio/public/phase1-data/ output/datasets/ --profile op2 --delete --dryrun

starter_dev <- c("s1d1", "s1d2")
train <- c(
  "s1d1", "s1d3",
  "s2d1", "s2d4", "s2d5",
  "s3d3", "s3d6", "s3d10"
)
valid <- c("s1d2", "s3d7")
test <- c("s4d1", "s4d8", "s4d9")

resources_test <- "resources_test/common/"
output_dir_p1v2 <- "output/datasets_2021-11-08/phase1v2/common/"
output_dir_p2 <- "output/datasets_2021-11-08/phase2_private/common/"


## process previous samplings
# prev <- anndata::read_h5ad("output/datasets/predict_modality/openproblems_bmmc_multiome_phase1_rna/openproblems_bmmc_multiome_phase1_rna.censor_dataset.output_train_mod2.h5ad", backed = "r")
# readr::write_lines(prev$var_names, "src/common/datasets/process_inhouse_datasets/sample_pm_atac_varnames.txt")

#############################
# Process multiome

# read and subset mod1 data
ad1old <- read_h5ad("output/manual_formatting/multiome/Multiome_gex_processed_train.h5ad", backed = "r")
ad1tr <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_gex_processed_fullfeat_train.h5ad")
ad1te <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_gex_processed_fullfeat_test_donors.h5ad")
ad1 <- anndata::concat(list(ad1tr[, colnames(ad1old)], ad1te[, colnames(ad1old)]))
ad1$obs <- ad1$obs %>% mutate_if(is.character, factor)
ad1$var <- ad1tr$var[colnames(ad1old), c("gene_ids-s1d1", "feature_types-s1d1", "genome-s1d1"), drop = FALSE]
colnames(ad1$var) <- c("gene_ids", "feature_types", "genome")
ad1$uns <- ad1tr$uns
rm(ad1tr, ad1te)

# read and subset mod2 data
ad2old <- read_h5ad("output/manual_formatting/multiome/Multiome_atac_processed_train.h5ad", backed = "r")
ad2tr <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_atac_processed_fullfeat_train.h5ad")
ad2te <- read_h5ad("output/manual_formatting_2021-11-08/atac/Multiome_atac_processed_fullfeat_test_donors.h5ad")
ad2 <- anndata::concat(list(ad2tr[, colnames(ad2old)], ad2te[, colnames(ad2old)]))
ad2$obs <- ad2$obs %>% mutate_if(is.character, factor)
ad2$uns <- ad2tr$uns
ad2$var <- ad2tr$var[colnames(ad2old), "feature_types-s1d1", drop = FALSE]
colnames(ad2$var) <- c("feature_types")
rm(ad2tr, ad2te)

adid <- "openproblems_bmmc_multiome"
assert_that(all(rownames(ad1) == rownames(ad2)))

unique(ad1$obs$batch)
assert_that(length(setdiff(ad1$obs$batch, c(train, valid, test))) == 0)
sort(setdiff(c(train, valid, test), unique(ad1$obs$batch)))

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
ad2$obsm$gene_activity <- as(ad2$obsm$gene_activity, "CsparseMatrix")

assert_that(max(ad1$X) < 100)
assert_that(max(ad2$X) < 100)

# copy over data
ad1$obs[["pseudotime_order_ATAC"]] <- ad2$obs[["pseudotime_order_ATAC"]]
ad2$obs[["pseudotime_order_GEX"]] <- ad1$obs[["pseudotime_order_GEX"]]

# # create samplings
ad2$uns$sample_pm_varnames <- readr::read_lines("src/common/datasets/process_inhouse_datasets/sample_pm_atac_varnames.txt")

# ad1$uns$sample_pm_obs # sample 1000
# ad2$uns$sample_pm_obs # sample 1000
# ad1$uns$sample_mm_train_obs # shuff
# ad1$uns$sample_mm_test_obs # shuff

#############
# create phase 1

ad1_phase1 <- ad1[ad1$obs$batch %in% c(train, valid), ]$copy()
ad2_phase1 <- ad2[ad2$obs$batch %in% c(train, valid), ]$copy()

ad1_phase1$obs[["is_train"]] <- ad1_phase1$obs[["batch"]] %in% train
ad2_phase1$obs[["is_train"]] <- ad1_phase1$obs[["batch"]] %in% train

adid_phase1 <- paste0(adid, "_phase1v2")
ad1_phase1$uns[["dataset_id"]] <- adid_phase1
ad2_phase1$uns[["dataset_id"]] <- adid_phase1

dir.create(paste0(output_dir_p1v2, adid_phase1), recursive = TRUE)
ad1_phase1$write_h5ad(paste0(output_dir_p1v2, adid_phase1, "/", adid_phase1, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_phase1$write_h5ad(paste0(output_dir_p1v2, adid_phase1, "/", adid_phase1, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")
rm(ad1_phase1, ad2_phase1)

#############
# create phase 2

ad1_phase2 <- ad1$copy()
ad2_phase2 <- ad2$copy()

ad1_phase2$obs[["is_train"]] <- ad1_phase2$obs[["batch"]] %in% c(train, valid)
ad2_phase2$obs[["is_train"]] <- ad1_phase2$obs[["batch"]] %in% c(train, valid)

adid_phase2 <- paste0(adid, "_phase2")
ad1_phase2$uns[["dataset_id"]] <- adid_phase2
ad2_phase2$uns[["dataset_id"]] <- adid_phase2

dir.create(paste0(output_dir_p2, adid_phase2), recursive = TRUE)
ad1_phase2$write_h5ad(paste0(output_dir_p2, adid_phase2, "/", adid_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
ad2_phase2$write_h5ad(paste0(output_dir_p2, adid_phase2, "/", adid_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")
rm(ad1_phase2, ad2_phase2)

#############
# create openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter

starter_cell <- which(ad1$obs$batch %in% starter_dev)
row_sel <- sort(sample(starter_cell, 500))
col_sel1 <- seq_len(ncol(ad1)) %in% sample.int(ncol(ad1), 600, prob = Matrix::colSums(ad1$X > 0))
col_sel2 <- seq_len(ncol(ad2)) %in% sample.int(ncol(ad2), 600, prob = Matrix::colSums(ad2$X > 0))

ad1_starter <- ad1[row_sel, col_sel1]$copy()
ad2_starter <- ad2[row_sel, col_sel2]$copy()

ad1_starter$obsm <- NULL
ad2_starter$obsm <- list(gene_activity = ad2_starter$obsm[["gene_activity"]])

ad1_starter$obs[["is_train"]] <- ad1_starter$obs[["batch"]] %in% train
ad2_starter$obs[["is_train"]] <- ad2_starter$obs[["batch"]] %in% train

adid_starter <- paste0(adid, "_starter")
ad1_starter$uns[["dataset_id"]] <- adid_starter
ad2_starter$uns[["dataset_id"]] <- adid_starter

dir.create(paste0(resources_test, adid_starter), recursive = TRUE)
ad1_starter$write_h5ad(paste0(resources_test, adid_starter, "/", adid_starter, ".output_rna.h5ad"), compression = "gzip")
ad2_starter$write_h5ad(paste0(resources_test, adid_starter, "/", adid_starter, ".output_mod2.h5ad"), compression = "gzip")
rm(ad1_starter, ad2_starter)



#############################
# Process cite

# read and subset mod1 data
bd1old <- read_h5ad("output/manual_formatting/cite/Cite_gex_processed_train.h5ad", backed = "r")
bd1tr <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_gex_processed_train_fullfeat.h5ad")
bd1te <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_gex_processed_test_donors_fullfeat.h5ad")
bd1 <- anndata::concat(list(bd1tr[, colnames(bd1old)], bd1te[, colnames(bd1old)]))
bd1$obs <- bd1$obs %>% mutate_if(is.character, factor)
bd1$var <- bd1tr$var[colnames(bd1old), , drop = FALSE] %>% select(-contains("-s4d"))
bd1$uns <- bd1tr$uns
rm(bd1tr, bd1te)

# read and subset mod2 data
bd2old <- read_h5ad("output/manual_formatting/cite/Cite_adt_processed_train.h5ad", backed = "r")
bd2tr <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_adt_processed_train_fullfeat.h5ad")
bd2te <- read_h5ad("output/manual_formatting_2021-11-08/cite/Cite_adt_processed_test_donors_fullfeat.h5ad")
bd2 <- anndata::concat(list(bd2tr[, colnames(bd2old)], bd2te[, colnames(bd2old)]))
bd2$obs <- bd2$obs %>% mutate_if(is.character, factor)
bd2$uns <- bd2tr$uns
bd2$var <- bd2tr$var[colnames(bd2old), , drop = FALSE] %>% select(-contains("-s4d"))
rm(bd2tr, bd2te)


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

bdid_phase1 <- paste0(bdid, "_phase1v2")
bd1_phase1$uns[["dataset_id"]] <- bdid_phase1
bd2_phase1$uns[["dataset_id"]] <- bdid_phase1

dir.create(paste0(output_dir_p1v2, bdid_phase1), recursive = TRUE)
bd1_phase1$write_h5ad(paste0(output_dir_p1v2, bdid_phase1, "/", bdid_phase1, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
bd2_phase1$write_h5ad(paste0(output_dir_p1v2, bdid_phase1, "/", bdid_phase1, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")
rm(bd1_phase1, bd2_phase1)

#############
# create phase 2

bd1_phase2 <- bd1$copy()
bd2_phase2 <- bd2$copy()

bd1_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)
bd2_phase2$obs[["is_train"]] <- bd1_phase2$obs[["batch"]] %in% c(train, valid)

bdid_phase2 <- paste0(bdid, "_phase2")
bd1_phase2$uns[["dataset_id"]] <- bdid_phase2
bd2_phase2$uns[["dataset_id"]] <- bdid_phase2

dir.create(paste0(output_dir_p2, bdid_phase2), recursive = TRUE)
bd1_phase2$write_h5ad(paste0(output_dir_p2, bdid_phase2, "/", bdid_phase2, ".manual_formatting.output_rna.h5ad"), compression = "gzip")
bd2_phase2$write_h5ad(paste0(output_dir_p2, bdid_phase2, "/", bdid_phase2, ".manual_formatting.output_mod2.h5ad"), compression = "gzip")
rm(bd1_phase2, bd2_phase2)


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
rm(bd1_starter, bd2_starter)

# compare sizes

orig <- c(
  list.files("output/datasets/common", recursive = TRUE, full.names = TRUE),
  list.files("output/datasets_2021-11-08/phase1v2/common", recursive = TRUE, full.names = TRUE),
  list.files("output/datasets_2021-11-08/phase2_private/common", recursive = TRUE, full.names = TRUE)
)

df <- map_df(
  orig,
  function(fn) {
    ad <- anndata::read_h5ad(fn, backed = "r")
    tibble(
      path = fn,
      dataset_id = ad$uns[["dataset_id"]],
      n_obs = nrow(ad),
      n_vars = ncol(ad),
      modality = unique(ad$var$feature_types),
      organism = ad$uns[["organism"]]
    )
  }
)

df %>%
  mutate(
    phase = gsub(".*phase([^_]*)", "\\1", dataset_id),
    platform = gsub(".*bmmc_([^_]*)_.*", "\\1", dataset_id)
  ) %>%
  arrange(platform, modality, phase) %>%
  select(-path)

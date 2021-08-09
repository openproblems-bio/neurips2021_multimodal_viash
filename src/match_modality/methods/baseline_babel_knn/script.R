cat(">> Loading dependencies\n")

requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("DropletUtils", quietly = TRUE)

options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input_mod1 = "output/task1_datasets/10x_pbmc_granulocyte_sorted_3k_rna/10x_pbmc_granulocyte_sorted_3k_rna.prepare_task1_dataset.output_mod1.h5ad",
  input_mod2 = "output/task1_datasets/10x_pbmc_granulocyte_sorted_3k_rna/10x_pbmc_granulocyte_sorted_3k_rna.prepare_task1_dataset.output_mod2.h5ad",
  output = "output_babel_knn.h5ad",
  data_id = "granulocyte_multiome10x"
  n_dims = 10,
  n_neigh= 10
)
if (!dir.exists(dirname(par$input_mod1))) {
  stop(
    "Directory '", dirname(par$input_mod1), "' does not exist.\n",
    "You probably need to run `aws s3 sync s3://neurips2021-multimodal-public-datasets/ output/` first."
  )
}


babel_location <- "../babel/bin/"      # location of babel executables

cat(">> Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
if (is.null(ad1$var$gene_ids)) ad1$var$gene_ids <- colnames(ad1)
ad2 <- anndata::read_h5ad(par$input_mod2)
if (is.null(ad2$var$gene_ids)) ad2$var$gene_ids <- colnames(ad2)

mod1 <- unique(ad1$var$feature_types)

# subset train and test
ad1_train <- ad1[ad1$obs$group == "train", ]
ad1_test <-  ad1[ad1$obs$group == "test", ]
# obs2 also has train/test split
ad2_train <- ad2[ad2$obs$group == "train", ]
ad2_test <-  ad2[ad2$obs$group == "test", ]

# multiome_matrix for export to Babel's input format
multiome_matrix <- cbind(ad1_train$X, ad2_train$X)



# generate multiome anndata object for training
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = bind_rows(ad1_train$var, ad2_train$var),
  obs = ad1_train$obs
)

# setting up babel dirs
tmpdir <- tempfile(pattern = "babel_temp", fileext = "/")
cat(">> Setting up directories for babel at ", tmpdir, "\n", sep = "")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE))

dir_data <- paste0(tmpdir, "data/")     # location of input files
dir.create(dir_data)
dir_model <- paste0(tmpdir, "model/")   # location of babel model
dir_pred <- paste0(tmpdir, "pred/")     # location of predictions

feature_type_map <- c(
  "GEX" = "Gene Expression",
  "ADT" = "Peaks", # try to make it run on ADT data as well
  "ATAC" = "Peaks"
)

cat(">> Writing train dataset as 10x-CellRanger H5 format\n")
DropletUtils::write10xCounts(
  paste0(dir_data, "train_input.h5"),
  t(ad_babel$X),
  gene.id = ad_babel$var$gene_ids,
  gene.symbol = colnames(ad_babel),
  barcodes = rownames(ad_babel),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[ad_babel$var$feature_types]),
  overwrite = TRUE
)

cat(">> Writing test dataset as 10x-CellRanger H5 format\n")
DropletUtils::write10xCounts(
  paste0(dir_data, "test_input.h5"),
  t(ad1_test$X),
  gene.id = ad1_test$var$gene_ids,
  gene.symbol = colnames(ad1_test),
  barcodes = rownames(ad1_test),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[ad1_test$var$feature_types]),
  overwrite = TRUE
)



cat(">> Babel: train model\n")
babel_train_cmd <- paste0(
  "/opt/conda/bin/conda run -n babel ",
  "python ", babel_location, "train_model.py ",
  "--data ", dir_data, "train_input.h5 ",
  "--outdir ", dir_model
)
out1 <- system(babel_train_cmd)

# check whether training succeeded
expect_equal(out1, 0, info = paste0("Model training failed with exit code ", out1))

cat(">> Babel: predict from model\n")
babel_pred_cmd <- paste0(
  "/opt/conda/bin/conda run -n babel ",
  "python ", babel_location, "predict_model.py ",
  "--checkpoint ", dir_model, " ",
  "--data ", dir_data, "test_input.h5 ",
  "--outdir ", dir_pred
)
out2 <- system(babel_pred_cmd)

# check whether prediction succeeded
expect_equal(out2, 0, info = paste0("Prediction failed with exit code ", out1))

cat(">> Read predictions\n")
pred <- anndata::read_h5ad(paste0(dir_pred, "/rna_atac_adata.h5ad"))


#######################################
#####  KNN

# Babel's output has filtered the features: select those in common
ad2_test_filter = ad2_test[,row.names(ad2_test$var) %in% row.names(pred$var)]

# train matrix for KNN
train_X = rbind(pred[,row.names(ad2_test_filter$var)]$X,ad2_test_filter$X)

# dimensional reduction
dr <- lmds::lmds(
  train_X,
  ndim = par$n_dims,
  distance_method = par$distance_method
)

# partition reduced matrix into train/test
n_cells = dim(dr)[1]

dr_knn = dr[ 1:n_cells , ] # predicted mod2 profiles babel
dr_y = dr[(n_cells+1):dim(dr)[1],] # actual mod2 profiles (shuffled)

# find k nearest neighbors for test data:
knn_indexes  = apply(dr_y, 1, function(x){
	FNN::get.knnx(dr_knn, query = x %>% t, k=par$n_neigh)$nn.index })

# matrix of matched cells
matched_cells = matrix(0, dim(knn_indexes)[2], dim(knn_indexes)[2])
# assing probabilities
for(i in 1:dim(matched_cells)[1])
	matched_cells[i, knn_indexes[,i]] = 1/par$n_neigh # assign equal probability to the k nearest neighbors

# make sparse matrix
out <- anndata::AnnData(
  X = as(matched_cells, "CsparseMatrix"),
  uns = list(
    dataset_id = ad1$uns[["dataset_id"]],
    method_id = "babel_knn"
  )
)

cat("Writing predictions to file\n")
zzz <- out$write_h5ad(par$output, compression = "gzip")

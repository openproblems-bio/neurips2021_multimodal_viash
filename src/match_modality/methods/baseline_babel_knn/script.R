cat(">> Loading dependencies\n")

requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(testthat, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("DropletUtils", quietly = TRUE)

options(tidyverse.quiet = TRUE)
library(tidyverse)

## VIASH START
par <- list(
  input_train_mod1 = "resources_test/match_modality/test_resource.train_mod1.h5ad",
  input_train_mod2 = "resources_test/match_modality/test_resource.train_mod2.h5ad",
  input_train_sol = "resources_test/match_modality/test_resource.train_sol.h5ad",
  input_test_mod1 = "resources_test/match_modality/test_resource.test_mod1.h5ad",
  input_test_mod2 = "resources_test/match_modality/test_resource.test_mod2.h5ad",
  output = "output.h5ad",
  n_dims = 10,
  n_neighs = 10
)
## VIASH END

cat("Reading h5ad files\n")
input_train_mod1 <- anndata::read_h5ad(par$input_train_mod1)
input_train_mod2 <- anndata::read_h5ad(par$input_train_mod2)
input_train_sol <- anndata::read_h5ad(par$input_train_sol)
input_test_mod1 <- anndata::read_h5ad(par$input_test_mod1)
input_test_mod2 <- anndata::read_h5ad(par$input_test_mod2)

babel_location <- "../babel/bin/"      # location of babel executables

cat(">> Reading h5ad files\n")
if (is.null(input_train_mod1$var$gene_ids)) input_train_mod1$var$gene_ids <- colnames(input_train_mod1)
if (is.null(input_train_mod2$var$gene_ids)) input_train_mod2$var$gene_ids <- colnames(input_train_mod2)
if (is.null(input_test_mod1$var$gene_ids)) input_test_mod1$var$gene_ids <- colnames(input_test_mod1)
if (is.null(input_test_mod2$var$gene_ids)) input_test_mod2$var$gene_ids <- colnames(input_test_mod2)

mod1 <- unique(input_test_mod1$var$feature_types)

# multiome_matrix for export to Babel's input format
multiome_matrix <- cbind(input_train_mod1$X, input_train_mod2$X)

# generate multiome anndata object for training
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = bind_rows(input_train_mod1$var, input_train_mod2$var)
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
  t(input_test_mod1$X),
  gene.id = input_test_mod1$var$gene_ids,
  gene.symbol = colnames(input_test_mod1),
  barcodes = rownames(input_test_mod1),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = unname(feature_type_map[input_test_mod1$var$feature_types]),
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
input_test_mod2_filter = input_test_mod2[,row.names(input_test_mod2$var) %in% row.names(pred$var)]

# train matrix for KNN
train_X = rbind(pred[,row.names(input_test_mod2_filter$var)]$X,input_test_mod2_filter$X)

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

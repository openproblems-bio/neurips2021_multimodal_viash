cat("Loading dependencies\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
requireNamespace("DropletUtils", quietly = TRUE)
# for "https://rdrr.io/github/barkasn/nbHelpers/man/transpose_dgRMatrix.html"
requireNamespace("nbHelpers", quietly = TRUE)

# TODO: Add setup for babel
# git clone https://github.com/wukevin/babel.git
# conda env create -f environment.yml

## VIASH START
par <- list(
  input_mod1 = "resources_test/task1/test_resource.mod1.h5ad",
  input_mod2 = "resources_test/task1/test_resource.mod2.h5ad",
  output = "output_babel.h5ad",
  data_id = "granulocyte_multiome10x"
)
## VIASH END

babel_location <- "../babel/bin/"     # location of babel executables

# folder names for Babel
tmpdir <- tempfile(pattern = "babel_temp")
dir.create(tmpdir)
on.exit(unlink(tmpdir, recursive = TRUE))

dir_data <- paste0(tmpdir, "/data")     # location of input files
dir.create(dir_data)
dir_model <- paste0(tmpdir, "/model")   # location of babel model
dir_pred <- paste0(tmpdir, "/pred") # location of predictions

cat("Reading h5ad files\n")
ad1 <- anndata::read_h5ad(par$input_mod1)
ad2 <- anndata::read_h5ad(par$input_mod2)

mod1 <- unique(ad1$var$feature_types)

# change feature types for each adata
if (mod1 == "GEX") {
  ad1$var$feature_types <- "Gene Expression"
  ad2$var$feature_types <- "Peaks"
} else {
  ad1$var$feature_types <- "Peaks"
  ad2$var$feature_types <- "Gene Expression"
}

# subset train and merge both categories (required by Babel)
ad1_train <- ad1[ad1$obs$group == "train", ]

# multiome_matrix for export to Babel's input format
multiome_matrix <- cbind(ad1_train$X, ad2$X)

# generate multiome anndata object
ad_babel <- anndata::AnnData(
  X = multiome_matrix,
  var = bind_rows(ad1_train$var, ad2$var),
  obs = ad1_train$obs
)

# write adata objects as 10x-CellRanger H5 format
# train object: test + train splits
DropletUtils::write10xCounts(
  paste0(dir_data, "/train_input.h5"),
  t(ad_babel$X),
  gene.id = row.names(ad_babel$var),
  gene.symbol = row.names(ad_babel$var),
  barcodes = row.names(ad_babel$obs),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = ad_babel$var$feature_types
)

# test object:
ad1_test <- ad1[ad1$obs$split == "test", ]

DropletUtils::write10xCounts(
  paste0(dir_data, "/test_input.h5"),
  t(ad1_test$X),
  gene.id = row.names(ad1_test$var),
  gene.symbol = row.names(ad1_test$var),
  barcodes = row.names(ad1_test$obs),
  type = "HDF5",
  version = "3",
  genome = "GRCh38",
  gene.type = ad1_test$var$feature_types
)

# RUN Babel
# switch conda environment as babel requires anndata v0.6
# TODO: set up conda environment in viash config.
use_condaenv("babel")

# train
babel_train_cmd <- paste0(
  "python ", babel_location, "/train_model.py ",
  "--data ", dir_data, "/",
  "train_input.h5 ",
  "--outdir ", dir_model
)
system(babel_train_cmd)
# TODO: check whether process failed or not

# test

babel_pred_cmd <- paste0(
  "python ", babel_location, "/predict_model.py ",
  "--checkpoint ", dir_model, " ",
  "--data ", dir_data, "/",
  "test_input.h5 ",
  "--outdir ", dir_pred
)
system(babel_pred_cmd)
# TODO: check whether process failed or not

# Babel generated a model folder with h5 files for the training set
# After prediction of the test data, the final h5 file is in ./babel_model_output
# format [input:RNA, prediction: ATAC]
out <- anndata::read_h5ad(paste0(dir_pred, "/rna_atac_adata.h5ad"))

# final h5 for prediction of test data
# add meta.data for the output
out$uns <- list(method_id = "babel")

# write output file
zzz <- out$write_h5ad(par$output, compression = "gzip")
